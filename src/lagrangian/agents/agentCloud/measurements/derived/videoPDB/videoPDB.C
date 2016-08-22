/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "videoPDB.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(videoPDB, 0);
addToRunTimeSelectionTable(agentMeasurement, videoPDB, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
videoPDB::videoPDB
(
    Time& t,
    const polyMesh& mesh,
    agentCloud& cloud,
    const dictionary& dict
)
:
    agentMeasurement(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    agentIds_(),
    fieldName_(propsDict_.lookup("fieldName")),
    iteration_(0),
    regionName_(),
    regionId_(-1),
    timeIndex_(0),    
    nSteps_(readLabel(propsDict_.lookup("numberOfOutputSteps"))),
    variableMols_(false),
    nSiteEstimate_(-1),
    startTime_(0.0),
    endTime_(GREAT),
    accumulatedTime_(0.0),
    deltaT_(t.deltaT().value())
{
    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();

    option_ = "mesh";
    
    if(propsDict_.found("option"))
    {
        const word option = propsDict_.lookup("option");
        option_ = option;
    }
    
    if(option_ == "zone")
    {
        const word regionName = propsDict_.lookup("zoneName");
        regionName_ = regionName;

        const cellZoneMesh& cellZones = mesh_.cellZones();

        regionId_ = cellZones.findZoneID(regionName_);

        if(regionId_ == -1)
        {
            FatalErrorIn("videoPDB::videoPDB()")
                << "Cannot find region: " << regionName_ << nl << "in: "
                << time_.time().system()/"fieldPropertiesDict"
                << exit(FatalError);
        }
    }
    
    if(option_ == "boundBox")
    {
        PtrList<entry> boxList(propsDict_.lookup("boxes"));

        boxes_.setSize(boxList.size());

        forAll(boxList, b)
        {
            const entry& boxI = boxList[b];
            const dictionary& dict = boxI.dict();

            vector startPoint = dict.lookup("startPoint");
            vector endPoint = dict.lookup("endPoint");
            boxes_[b].resetBoundedBox(startPoint, endPoint);
        }
    }
    
    scalingFactor_ = 1;
    
    if (propsDict_.found("scalingFactor"))
    {
        scalingFactor_ = readScalar(propsDict_.lookup("scalingFactor"));
    }   
//     if (propsDict_.found("molOption"))
//     {
//         const word molOption = propsDict_.lookup("molOption");
//         
//         molOption_ = molOption;
//     }
    
    if (propsDict_.found("variableMols"))
    {
        variableMols_ = Switch(propsDict_.lookup("variableMols"));
        
        nSiteEstimate_ = readLabel(propsDict_.lookup("nSiteEstimate"));
        rDummy_ = propsDict_.lookup("outsidePosition");
    }


    if (propsDict_.found("startAtTime"))
    {    
        startTime_ = readScalar(propsDict_.lookup("startAtTime"));
    }
    
    if (propsDict_.found("endAtTime"))
    {
        endTime_ = readScalar(propsDict_.lookup("endAtTime"));
    }
    
    writeFirstTimeStep_ = true;
    
    if (propsDict_.found("writeFirstTimeStep"))
    {    
        writeFirstTimeStep_ = readScalar(propsDict_.lookup("writeFirstTimeStep"));
    }    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

videoPDB::~videoPDB()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void videoPDB::createField()
{
    // set many files

    label nAgents = 0;
    
    {    
        IDLList<agent>::iterator mol(cloud_.begin());

        for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
        {
            if(findIndex(agentIds_, mol().id()) != -1)
            {
                nAgents++;
            }
        }    
    }
    
    if (Pstream::parRun())
    {
        reduce(nAgents, sumOp<label>());
    }
    
    n_ = label(nAgents/100000) + 1;

    if(n_ == 0) 
    {
        FatalErrorIn("videoPDB::videoPDB()")
            << " number of files should be at least 1." << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }
    else if (n_ == 1)
    {
        Info << "videoPDB" << nl
             << "-> number of files set to = " << n_ 
             << nl << endl;        
    }
    else
    {
        Info << "WARNING in videoPDB" << nl
             << "-> number of files set to = " << n_ 
             << nl << endl;
    }
    
    minLimit_.setSize(n_, -1);
    maxLimit_.setSize(n_, -1);
    minLimit_[0] = 0;
    maxLimit_[0] = 99999;

    for (int i = 1; i < n_; i++)
    {
        minLimit_[i] = 100000*(i);
        maxLimit_[i] = (100000*(i+1)) - 1;
    }
    
    
    if(writeFirstTimeStep_)
    {
        iteration_++;        
        write();
    }
}

void videoPDB::calculateField()
{
    
    accumulatedTime_ += deltaT_;
    
    if((accumulatedTime_ >= startTime_) && (accumulatedTime_ <= endTime_))
    {     
        timeIndex_++;
     
        if(timeIndex_ >= nSteps_)
        {
            iteration_++;

            write();

            timeIndex_ = 0;
        }
    }
}

void videoPDB::writeField()
{}

void videoPDB::writeInMesh(List<labelField>& agentIds, List<vectorField>& sites)
{
    label myProc =  Pstream::myProcNo();

    IDLList<agent>::iterator mol(cloud_.begin());

    DynamicList<vector> sitePositions(0);
    DynamicList<label> agentIDs(0);

    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {
        if(findIndex(agentIds_, mol().id()) != -1)
        {
            agentIDs.append(mol().id());
            sitePositions.append(mol().position());
        }
    }

    sites[myProc].transfer(sitePositions);
    agentIds[myProc].transfer(agentIDs);
}

void videoPDB::writeInBoundBox(List<labelField>& agentIds, List<vectorField>& sites)
{
    label myProc =  Pstream::myProcNo();

    IDLList<agent>::iterator mol(cloud_.begin());

    DynamicList<vector> sitePositions(0);
    DynamicList<label> agentIDs(0);

    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {
        forAll(boxes_, b)
        {
            if(boxes_[b].contains(mol().position()))
            {        
                if(findIndex(agentIds_, mol().id()) != -1)
                {
                    agentIDs.append(mol().id());
                    sitePositions.append(mol().position());
                }
            }
        }
    }

    sites[myProc].transfer(sitePositions);
    agentIds[myProc].transfer(agentIDs);
}

void videoPDB::writeInZone(List<labelField>& agentIds, List<vectorField>& sites)
{
    label myProc =  Pstream::myProcNo();

    {
        const List< DynamicList<agent*> >& cellOccupancy
            = cloud_.cellOccupancy();

        const labelList& cells = mesh_.cellZones()[regionId_];

        DynamicList<vector> sitePositions(0);
        DynamicList<label> agentIDs(0);

        forAll(cells, c)
        {
            const label& cellI = cells[c];
            const List<agent*>& molsInCell = cellOccupancy[cellI];

            forAll(molsInCell, mIC)
            {
                agent* molI = molsInCell[mIC];

                if(findIndex(agentIds_, molI->id()) != -1)
                {
                    agentIDs.append(molI->id());
                    sitePositions.append(molI->position());
                }
            }
        }

        sites[myProc].transfer(sitePositions);
        agentIds[myProc].transfer(agentIDs);
    }
}

void videoPDB::write()
{
    List<labelField> agentIds(Pstream::nProcs());
    List<vectorField> sites(Pstream::nProcs());

    label myProc =  Pstream::myProcNo();

    if(option_ == "zone")
    {
        Info << "videoPDB: write in zone" << endl;

        writeInZone(agentIds, sites);
    }
    if(option_ == "boundBox")
    {
        Info << "videoPDB: write in mesh" << endl;

        writeInBoundBox(agentIds, sites);       
    }
    if(option_ == "mesh")
    {
        Info << "videoPDB: write in mesh" << endl;

        writeInMesh(agentIds, sites);
    }

    label totalSites = sites[myProc].size();
    label totalMols = agentIds[myProc].size();

    if (Pstream::parRun())
    {
        reduce(totalSites, sumOp<label>());
        reduce(totalMols, sumOp<label>());
    }

    if (Pstream::parRun())
    {
        // send to master (master does not send)
        if(!Pstream::master())
        {
            const int proc = 0;
            {
                OPstream toNeighbour(Pstream::blocking, proc);
                toNeighbour << sites[myProc] << agentIds[myProc];
            }
        }

        //- receiving (master only receives)
        if(Pstream::master())
        {
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    vectorField sitesProc;
                    labelField agentIdsProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> sitesProc >> agentIdsProc;
                    }

                    sites[p].setSize(sitesProc.size());
                    agentIds[p].setSize(agentIdsProc.size());

                    forAll(sitesProc, i)
                    {
                        sites[p][i] = sitesProc[i];
                    }

                    forAll(agentIdsProc, i)
                    {
                        agentIds[p][i] = agentIdsProc[i];
                    }
                }
            }
        }
    }


    Info << "totalSites: " << totalSites << endl;

    scalar nFiles = totalSites/99999;

    if(nFiles > n_)
    {
        Info<< "WARNING. Error in PDB. You need a total number of files of: "
            << nFiles << " not: " << n_
            << endl;
    }

    if(Pstream::master())
    {
        const reducedUnits& rU = cloud_.redUnits();

        for (int j = 0; j < n_; j++)
        {
            std::string s;
            std::stringstream out;
            out << j;
            s = out.str();

            fileName fName(casePath_/"agentCloud_"+fieldName_+"_"+s+".pdb");
    
            std::ofstream file(fName.c_str(),ios_base::app);
        
            if(file.is_open())
            {
                file << "MODEL " << iteration_ << nl;
    
                label nSites = 0;
                label nMols = 1;
    
   
                // for all processors
                forAll(agentIds, p)
                {
                    label posCounter = -1;
                    
                    forAll(agentIds[p], i)
                    {
                        label id = agentIds[p][i];

                        nSites++;
                        posCounter++;
                        
                        if((nSites >= minLimit_[j] ) && (nSites <= maxLimit_[j]))
                        {
                            vector rS = sites[p][posCounter]*rU.refLength()*scalingFactor_;

                            // site 
                            file.width(6);
                            file << std::left << "ATOM";
                            file.width(5);
                            file << std::right << nSites-minLimit_[j];
                            file << "  ";
                            file.width(3);
                            file << std::left << cloud_.cP().agentIds()[id];
                            file << " ";
                            file.width(3);
                            file << std::right << cloud_.cP().agentIds()[id];
                            file << " ";
                            file.width(5);
                            file << nMols;
                            file << "    ";
                            file.width(8);
                            file.precision(3);
                            file.setf(std::ios::fixed,std::ios::floatfield);  
                            file << rS.x();
                            file.width(8);
                            file.precision(3);
                            file.setf(std::ios::fixed,std::ios::floatfield);  
                            file << rS.y();
                            file.width(8);
                            file.precision(3);
                            file.setf(std::ios::fixed,std::ios::floatfield);  
                            file << rS.z();
                            file << "  1.00  0.00 ";
                            file << nl;
                            
                            nMols++;
                        }

//                         if((nSites >= minLimit_[j] ) && (nSites <= maxLimit_[j]))
//                         {
//                             nMols++;
//                         }
                    }
                }
                
                if(variableMols_)
                {
                    
                    label nBufferMols = nSiteEstimate_ - nSites;
                    
                    if(nBufferMols < 0)
                    {
                        FatalErrorIn("void atomisticPDB::writeField()")
                            << "Exceeded limits of estimated nMol. Increase -> " << nSiteEstimate_
                            << abort(FatalError);
                    }               
                    
                    label id = agentIds_[0];
                    
                    for (int i = 0; i < nBufferMols; i++)
                    {
                        nSites++;

                        if((nSites >= minLimit_[j] ) && (nSites <= maxLimit_[j]))
                        {
                            vector rS = rDummy_*rU.refLength()*scalingFactor_;

                            file.width(6);
                            file << std::left << "ATOM";
                            file.width(5);
                            file << std::right << nSites-minLimit_[j];
                            file << "  ";
                            file.width(3);
                            file << std::left << cloud_.cP().agentIds()[id];
                            file << " ";
                            file.width(3);
                            file << std::right << cloud_.cP().agentIds()[id];
                            file << " ";
                            file.width(5);
                            file << nMols;
                            file << "    ";
                            file.width(8);
                            file.precision(3);
                            file.setf(std::ios::fixed,std::ios::floatfield);  
                            file << rS.x();
                            file.width(8);
                            file.precision(3);
                            file.setf(std::ios::fixed,std::ios::floatfield);  
                            file << rS.y();
                            file.width(8);
                            file.precision(3);
                            file.setf(std::ios::fixed,std::ios::floatfield);  
                            file << rS.z();
                            file << "  1.00  0.00 ";
                            file << nl;
                            
//                             if((nSites >= minLimit_[j] ) && (nSites <= maxLimit_[j]))
//                             {
                                nMols++;
//                             }
                        }
                    }
                }
                
                file << "ENDMDL" <<  nl;
            }
            else
            {
                FatalErrorIn("void combinedPDB::writeField()")
                    << "Cannot open file " << fName
                    << abort(FatalError);
            }
    
            file.close();
        }
    }
}

void videoPDB::measureDuringForceComputation
(
    agent* molI,
    agent* molJ
)
{}

void videoPDB::measureDuringForceComputationSite
(
    agent* molI,
    agent* molJ,
    label sI,
    label sJ
)
{}

// const propertyField& videoPDB::fields() const
// {
//     return fields_;
// }

} // End namespace Foam

// ************************************************************************* //
