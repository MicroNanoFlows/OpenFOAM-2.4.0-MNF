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

#include "dropletInfo.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dropletInfo, 0);

addToRunTimeSelectionTable(utilField, dropletInfo, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dropletInfo::dropletInfo
(
    Time& t,
    const polyMesh& mesh,
    const reducedUnits& rU,
    const dictionary& dict
)
:
    utilField(t, mesh, rU, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.lookup("fieldName")),
    nameFile1_(propsDict_.lookup("nameFile_positions")),
    nameFile2_(propsDict_.lookup("nameFile_energies")),
    nameFile3_(propsDict_.lookup("nameFile_forces")),
    nameFile4_(propsDict_.lookup("nameFile_nPairs")),
    nameFile5_(propsDict_.lookup("nameFile_referenceMols")),
    endTime_(readScalar(propsDict_.lookup("endTime"))),
    mass_(readScalar(propsDict_.lookup("molecularMass")))
{
//     writeTime_ = GET;
    deltaT_ = time_.deltaT().value();
    
//     endTime_ = time_.endTime().value();
    
//     Info << "endTime = " << endTime_ << endl;
    
    
    mass_ /= rU_.refMass();
    
//     nSteps_ = endTime_/deltaT_;
    
//     nOutputSteps_=1000;
//     
//     if(propsDict_.found("nOutputSteps"))
//     {
//         nOutputSteps_ = readLabel(propsDict_.lookup("nOutputSteps"));
//     }
//     
    startTime_ = deltaT_*2.0;
    
    if(propsDict_.found("startTime"))    
    {        
        startTime_ = readScalar(propsDict_.lookup("startTime"));
    }
//     writeTime/deltaT
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dropletInfo::~dropletInfo()
{}





// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dropletInfo::createField()
{
    
}


void dropletInfo::findSegments(List<DynamicList<label> >& ids)
{
    // outer list = number of segments
    // inner list = mol indices within segment 
    DynamicList<vector> segmentPositions;
    
    label counter = 0;
    
    forAll(positions_, i)
    {
        const vector& pI = positions_[i];
        
        label N = segmentPositions.size();
        
        if(fluidMols_[i])
        {
            if(N == 0)
            {
                segmentPositions.append(pI);
                ids[counter].append(i);
                counter++;
            }
            else
            {
                bool segmentFound = false;
                
                forAll(segmentPositions, j)
                {
                    const vector& pJ = segmentPositions[j];
                    
                    scalar pIJ = mag(pI - pJ);
                    
                    if(pIJ < D_)
                    {
                        segmentFound = true;
                        ids[j].append(i);
                    }
                }
                
                if(!segmentFound)
                {
                    segmentPositions.append(pI);
                    ids[counter].append(i);          
                    counter++;
                }
            }
        }
    }

    Info << "segmentPositions = " << segmentPositions << endl; 
}


void dropletInfo::getError
(
    const List<DynamicList<label> >& ids,
    label& error,
    label& id
)
{
    label actualN = nActualMols_;
    error = actualN;
    
    id = -1;
    
    forAll(ids, i)
    {
        label nMols = ids[i].size();
        
        if((actualN - nMols) < error)
        {
            error = actualN - nMols;
            id = i;
        }
    }
    
    Info << "error = " << error << endl;
}

void dropletInfo::setDropletPositions
(
    const List<DynamicList<label> >& ids,
    const label& id
)
{
    newMolsToPick_.setSize(fluidMols_.size(), false);
    
    nSelectedMols_ = 0;
    forAll(ids[id], j)
    {
        const label& tN = ids[id][j];
        newMolsToPick_[tN] = true;
        nSelectedMols_++;
    }
    
}


void dropletInfo::calculateField()
{
    if(time_.timeOutputValue() >= startTime_) // skipping some known issues in first time step
    {
        readFromStorage();
    
        {
            // estimate radius             
            label N = 0;
            
            forAll(positions_, i)
            {
                if(fluidMols_[i])
                {
                    N++;
                }
            }
            
            nActualMols_ = N;
            
            scalar rho = 1000.0/rU_.refMassDensity();
            scalar factor = 0.5;
            scalar term = (3.0*mass_*scalar(N))/(factor*rho*4.0*constant::mathematical::pi);
            scalar R = Foam::pow(term, (1.0/3.0));
            
            Info << " R = " << R << endl;
            
            scalar offset_ = 1.0;
            
            D_ = 2*R + offset_;
            
            boundedBox bbMesh( mesh_.bounds().min(), mesh_.bounds().max());
            scalar Lx = bbMesh.span().x();
            scalar Lz = bbMesh.span().z();
            scalar midX = Lx*0.5;            
            scalar midZ = Lz*0.5;            
            
            
            {

                List<DynamicList<label> > ids(10);
                
                findSegments(ids);
               
                bool fullSegmentFound = false;
                
                // try 0

                label error;
                label id;
                
                getError(ids, error, id);
                
                if(error < 5)
                {
                    fullSegmentFound = true;
                    
                    setDropletPositions(ids, id);
                }
            
            
                if(!fullSegmentFound)
                {
                    // repeat but first translate some molecules
                    forAll(ids, i)
                    {
                        //push all molecules in +Lx
                        // then push all molecules in +Lz
                        forAll(ids[i], j)
                        {
                            const label& tN = ids[i][j];
                        
                            vector r1 = positions_[tN];
                            
                            if(r1.x() < midX)
                            {
                                vector r2 = r1 + Lx*vector(1, 0, 0);
                                positions_[tN] = r2;
                            }
                        }
                    }
                    
                    // try again 
                    
                    List<DynamicList<label> > ids1(10);
            
                    findSegments(ids1);

                    label error1;
                    label id1;
                    
                    getError(ids1, error1, id1);
                    
                    if(error1 < 5)
                    {
                        fullSegmentFound = true;
                        setDropletPositions(ids1, id1);
                    }                        
                    
                    if(!fullSegmentFound)
                    {
                        forAll(ids1, i)
                        {
                            //push all molecules in +Lx
                            // then push all molecules in +Lz
                            forAll(ids1[i], j)
                            {
                                const label& tN = ids1[i][j];
                            
                                vector r1 = positions_[tN];
                                
                                if (r1.z() < midZ) 
                                {
                                    vector r2 = r1 + Lz*vector(0, 0, 1);
                                    positions_[tN] = r2;
                                }
                            }
                        }
                        
                        // try again
                        List<DynamicList<label> > ids2(10);
                
                        findSegments(ids2);

                        label error2;
                        label id2;
                        
                        getError(ids2, error2, id2);
                        
                        if(error2 < 5)
                        {
                            fullSegmentFound = true;
                            setDropletPositions(ids2, id2);
                        }        
                        else
                        {
                            FatalErrorIn("void dropletInfo::calculateField()")
                                << "Something strange - bad coding or something"
                                << exit(FatalError);                                
                        }
                    }
                }
            }
        }
    }
}
    
void dropletInfo::writeField()
{
    // write xmol file
    {
        const reducedUnits& rU = rU_;

        fileName fName(outputPath_/"shifted_"+fieldName_+".xmol");
    
        OFstream str(fName);
    
        str << nSelectedMols_ << nl << "polyMoleculeCloud site positions in angstroms" << nl;

        // for all processors
        forAll(positions_, i)
        {
            vector rI = positions_[i]*rU.refLength()*1e10;
            
            if(newMolsToPick_[i])
            {
                str <<  'C'
                    << ' ' << rI.x()
                    << ' ' << rI.y()
                    << ' ' << rI.z()
                    << nl;
            }
        }
    }    
    
    if(time_.timeOutputValue() >= endTime_)
    {
        
        
    }
/*    writeTimeData
    (
        outputPath_,
        "dropletInfo_"+fieldName_+"_method_1_cumulFreePath_nm_vs_time_ns.txt",
        time*rU_.refTime()/1e-9,
        mfp*rU_.refLength()/1e-9,
        false
    ); */    
}





void dropletInfo::readFromStorage()
{
    Info << "reading from storage" << endl;
    
    fileName pathName(inputPath_);
       
    {
        IFstream file(pathName/nameFile1_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> positions_;
        }
    }
    {
        IFstream file(pathName/nameFile2_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> energies_;
        }
    }    
    {
        IFstream file(pathName/nameFile3_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> forces_;
        }
    }
    {
        IFstream file(pathName/nameFile4_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> nPairs_;
        }
    }
    {
        IFstream file(pathName/nameFile5_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> fluidMols_;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
