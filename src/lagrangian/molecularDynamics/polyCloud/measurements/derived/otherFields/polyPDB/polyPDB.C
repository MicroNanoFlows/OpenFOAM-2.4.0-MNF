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

#include "polyPDB.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyPDB, 0);
addToRunTimeSelectionTable(polyField, polyPDB, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyPDB::polyPDB
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fields_(t, mesh, "dummy"),
    molIds_(),
    excludeSites_(),
    fieldName_(propsDict_.lookup("fieldName")),
    n_(readLabel(propsDict_.lookup("numberOfFiles"))),
    minLimit_(n_, -1),
    maxLimit_(n_, -1),
    iteration_(0),
    zone_(false),
    regionName_(),
    regionId_(-1),
    timeIndex_(0),    
    nSteps_(readLabel(propsDict_.lookup("numberOfOutputSteps"))),
    variableMols_(false),
    nMolsEstimate_(-1),
    startTime_(0.0),
    endTime_(GREAT),
    accumulatedTime_(0.0),
    deltaT_(t.deltaT().value())
{
    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    if(propsDict_.found("zoneName"))
    {
        zone_ = true;

        const word regionName = propsDict_.lookup("zoneName");
        regionName_ = regionName;

        const cellZoneMesh& cellZones = mesh_.cellZones();

        regionId_ = cellZones.findZoneID(regionName_);

        if(regionId_ == -1)
        {
            FatalErrorIn("polyPDB::polyPDB()")
                << "Cannot find region: " << regionName_ << nl << "in: "
                << time_.time().system()/"fieldPropertiesDict"
                << exit(FatalError);
        }
    }

    if(n_ == 0) 
    {
        FatalErrorIn("polyPDB::polyPDB()")
            << " number of files should be at least 1." << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }
    
    if (propsDict_.found("variableMols"))
    {
        variableMols_ = Switch(propsDict_.lookup("variableMols"));
        
        nMolsEstimate_ = readLabel(propsDict_.lookup("molEstimate"));
        rDummy_ = propsDict_.lookup("outsidePosition");
    }

    minLimit_[0] = 0;
    maxLimit_[0] = 99999;

    for (int i = 1; i < n_; i++)
    {
        minLimit_[i] = 100000*(i);
        maxLimit_[i] = (100000*(i+1)) - 1;
    }
    
    if (propsDict_.found("startAtTime"))
    {    
        startTime_ = readScalar(propsDict_.lookup("startAtTime"));
    }
    
    if (propsDict_.found("endAtTime"))
    {
        endTime_ = readScalar(propsDict_.lookup("endAtTime"));
    }    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyPDB::~polyPDB()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyPDB::createField()
{
    selectSiteIds sites
    (
        molCloud_.cP(),
        propsDict_,
        "sitesToExclude"
    );

    List<word> siteNames = sites.siteIdNames();

    excludeSites_.transfer(siteNames);

    Info   << "sites to exclude: " << excludeSites_ << endl;
}

void polyPDB::calculateField()
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

void polyPDB::writeField()
{}

void polyPDB::writeInMesh(List<labelField>& molIds, List<vectorField>& sites)
{
    label myProc =  Pstream::myProcNo();

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    DynamicList<vector> sitePositions(0);
    DynamicList<label> moleculeIds(0);

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            moleculeIds.append(mol().id());

//             const polyMolecule::constantProperties cP(molCloud_.constProps(mol().id()));

            forAll(mol().sitePositions(), i)
            {
                if(findIndex(excludeSites_, molCloud_.cP().siteNames(mol().id())[i]) == -1)
                {
                    sitePositions.append(mol().sitePositions()[i]);
                }
            }
        }
    }

    //sites[myProc].transfer(sitePositions.shrink());
    //molIds[myProc].transfer(moleculeIds.shrink());

    sites[myProc].transfer(sitePositions);
    molIds[myProc].transfer(moleculeIds);
}

void polyPDB::writeInZone(List<labelField>& molIds, List<vectorField>& sites)
{
    label myProc =  Pstream::myProcNo();

    {
        const List< DynamicList<polyMolecule*> >& cellOccupancy
            = molCloud_.cellOccupancy();

        const labelList& cells = mesh_.cellZones()[regionId_];

        DynamicList<vector> sitePositions(0);
        DynamicList<label> moleculeIds(0);

        forAll(cells, c)
        {
            const label& cellI = cells[c];
            const List<polyMolecule*>& molsInCell = cellOccupancy[cellI];

            forAll(molsInCell, mIC)
            {
                polyMolecule* molI = molsInCell[mIC];

                if(findIndex(molIds_, molI->id()) != -1)
                {
                    moleculeIds.append(molI->id());
//                     const polyMolecule::constantProperties cP(molCloud_.constProps(molI->id()));

                    forAll(molI->sitePositions(), i)
                    {
                        if(findIndex(excludeSites_, molCloud_.cP().siteNames(molI->id())[i]) == -1)
                        {
                            sitePositions.append(molI->sitePositions()[i]);
                        }
                    }
                }
            }
        }

        //sites[myProc].transfer(sitePositions.shrink());
        //molIds[myProc].transfer(moleculeIds.shrink());

        sites[myProc].transfer(sitePositions);
        molIds[myProc].transfer(moleculeIds);
    }
}

void polyPDB::write()
{
    List<labelField> molIds(Pstream::nProcs());
    List<vectorField> sites(Pstream::nProcs());

    label myProc =  Pstream::myProcNo();

    if(zone_)
    {
        Info << "polyPDB: write in zone" << endl;

        writeInZone(molIds, sites);
    }
    else
    {
        Info << "polyPDB: write in mesh" << endl;

        writeInMesh(molIds, sites);
    }

    label totalSites = sites[myProc].size();
    label totalMols = molIds[myProc].size();

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
                toNeighbour << sites[myProc] << molIds[myProc];
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
                    labelField molIdsProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> sitesProc >> molIdsProc;
                    }

                    sites[p].setSize(sitesProc.size());
                    molIds[p].setSize(molIdsProc.size());

                    forAll(sitesProc, i)
                    {
                        sites[p][i] = sitesProc[i];
                    }

                    forAll(molIdsProc, i)
                    {
                        molIds[p][i] = molIdsProc[i];
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
        const reducedUnits& rU = molCloud_.redUnits();

        for (int j = 0; j < n_; j++)
        {
            std::string s;
            std::stringstream out;
            out << j;
            s = out.str();

            fileName fName(casePath_/"polyMoleculeCloud_"+fieldName_+"_"+s+".pdb");
    
            std::ofstream file(fName.c_str(),ios_base::app);
        
            if(file.is_open())
            {
                file << "MODEL " << iteration_ << nl;
    
                label nSites = 0;
                label nMols = 1;
    
                const List<word>& idList(molCloud_.cP().molIds());
    
                // for all processors
                forAll(molIds, p)
                {
                    label posCounter = -1;
    
                    forAll(molIds[p], i)
                    {
                        label molId = molIds[p][i];
    
//                         const polyMolecule::constantProperties cP(molCloud_.constProps(molId));
    
                        label n = molCloud_.cP().nSites(molId);

                        for (int i = 0; i < n; i++)
                        {
                            if(findIndex(excludeSites_,  molCloud_.cP().siteNames(molId)[i]) == -1)
                            {
                                nSites++;
                                posCounter++;
    
                                if((nSites >= minLimit_[j] ) && (nSites <= maxLimit_[j]))
                                {
                                    vector rS = sites[p][posCounter]*rU.refLength()*1.0e10;

                                    if(molCloud_.cP().pointMolecule(molId))
                                    {
                                        // site H1
                                        file.width(6);
                                        file << std::left << "ATOM";
                                        file.width(5);
                                        file << std::right << nSites-minLimit_[j];
                                        file << "  ";
                                        file.width(3);
                                        file << std::left << molCloud_.cP().siteNames(molId)[i];
                                        file << " ";
                                        file.width(3);
                                        file << std::right << molCloud_.cP().siteNames(molId)[i];
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
                                    }
                                    else
                                    {
                                        if(idList[molId] == "water")
                                        {
                                            file << "HETATM";
                                            file.width(5);
                                            file << nSites-minLimit_[j];
                                            file << "  ";
                                            file.width(3);
                                            file << std::left << molCloud_.cP().siteNames(molId)[i];
                                            file << " ";
                                            file.width(3);
                                            file << std::right << "HOH";
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
                                        }
            
                                        else
                                        {
                                            file << "HETATM";
                                            file.width(5);
                                            file << nSites-minLimit_[j];
                                            file << "  ";
                                            file.width(3);
                                            file << std::left << molCloud_.cP().siteNames(molId)[i];
                                            file << " ";
                                            file.width(3);
                                            file << std::right << "XXX";
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
                                        }
                                    }
                                }
                            }
                        }

                        if((nSites >= minLimit_[j] ) && (nSites <= maxLimit_[j]))
                        {
                            nMols++;
                        }
                    }
                }
                
                if(variableMols_)
                {
                    label nBufferMols = nMolsEstimate_ - nSites;
                    
                    if(nBufferMols < 0)
                    {
                        FatalErrorIn("void combinedPDB::writeField()")
                            << "Exceeded limits of estimated nMol. Increase -> " << nMolsEstimate_
                            << ", to at least -> " << nSites 
                            << abort(FatalError);
                    }               
                    
                    label molId = molIds_[0];
//                     const polyMolecule::constantProperties cP(molCloud_.constProps(molId));
                    word siteName = molCloud_.cP().siteNames(molId)[0];
                    
                    for (int i = 0; i < nBufferMols; i++)
                    {
                        nSites++;

                        if((nSites >= minLimit_[j] ) && (nSites <= maxLimit_[j]))
                        {
                            vector rS = rDummy_*rU.refLength()*1.0e10;

                            if(molCloud_.cP().pointMolecule(molId))
                            {
                                // site H1
                                file.width(6);
                                file << std::left << "ATOM";
                                file.width(5);
                                file << std::right << nSites-minLimit_[j];
                                file << "  ";
                                file.width(3);
                                file << std::left << molCloud_.cP().siteNames(molId)[0];
                                file << " ";
                                file.width(3);
                                file << std::right << molCloud_.cP().siteNames(molId)[0];
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
                            }
                            else
                            {

                                file << "HETATM";
                                file.width(5);
                                file << nSites-minLimit_[j];
                                file << "  ";
                                file.width(3);
                                file << std::left << molCloud_.cP().siteNames(molId)[0];
                                file << " ";
                                file.width(3);
                                file << std::right << "XXX";
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

                            }
                            
                            nMols++;
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

void polyPDB::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyPDB::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& polyPDB::fields() const
{
    return fields_;
}

} // End namespace Foam

// ************************************************************************* //
