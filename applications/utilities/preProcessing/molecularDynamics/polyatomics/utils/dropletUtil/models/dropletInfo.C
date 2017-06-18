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
    translate_(propsDict_.lookup("translate")),
    endTime_(readScalar(propsDict_.lookup("endTime"))),
    maxR_(readScalar(propsDict_.lookup("maxR")))
{
//     writeTime_ = GET;
    deltaT_ = time_.deltaT().value();
    
//     endTime_ = time_.endTime().value();
    
//     Info << "endTime = " << endTime_ << endl;
    
    
//     mass_ /= rU_.refMass();
    
//     nSteps_ = endTime_/deltaT_;
    
//     nOutputSteps_=1000;
//     
//     if(propsDict_.found("nOutputSteps"))
//     {
//         nOutputSteps_ = readLabel(propsDict_.lookup("nOutputSteps"));
//     }
//     
    nSegs_ = 5;
    
    if(propsDict_.found("nSegments"))
    {
        nSegs_ = readLabel(propsDict_.lookup("nSegments"));
    }    
    
    startTime_ = deltaT_*2.0;
    
    if(propsDict_.found("startTime"))    
    {        
        startTime_ = readScalar(propsDict_.lookup("startTime"));
    }
    
    val_ = 1;
    
    boundedBox bbMesh( mesh_.bounds().min(), mesh_.bounds().max() );
    scalar Lx = bbMesh.span().x();
    scalar Ly = bbMesh.span().x();
    scalar Lz = bbMesh.span().z();

    
    if(translate_ == "X")
    {
        translationVector_ = vector(1, 0, 0)*Lx;
    }
    
    if(translate_ == "Y")
    {
        translationVector_ = vector(0, 1, 0)*Ly;
    }       
    
    if(translate_ == "Z")
    {
        translationVector_ = vector(0, 0, 1)*Lz;
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

void dropletInfo::reconstructDroplet()
{
    Info << " reconstruct droplet ... " << endl;
    
    {
        label N = 0;

        forAll(positions_, i)
        {
            if(fluidMols_[i])
            {
                N++;
            }
        }
        
        nFluidMols_ = N;

        Info << "no. of fluid mols = " << N << endl;
        
        List<vector> fluidPositions(N);
        List<label> fluidIds(N);
        
        label c = 0;
        
        forAll(positions_, i)
        {
            if(fluidMols_[i])
            {
                fluidPositions[c] = positions_[i];
                fluidIds[c] = i;
                c++;
            }
        }            
        
        List<DynamicList<vector> > selectPositions(nSegs_);
        List<DynamicList<label> > selectIds(nSegs_);

        // initial estimate
        
        List<label> molsChosen(N, -1);
        
        forAll(selectPositions, s)
        {
            vector com = vector::zero;
            
            if(s == 0)
            {
                com = fluidPositions[0];                
            }
            else
            {
                bool foundMol = false;
                
                forAll(molsChosen, i)
                {
                    if(molsChosen[i] == -1)
                    {
                        com = fluidPositions[i];
                        foundMol = true;
                        break;
                    }
                }
                
                if(!foundMol)
                {
                    break;
                }
            }
            
            scalar count = 1.0;
            
            label molsCount = 0;
            
            bool finish = false;
            
            while(!finish)
            {
                forAll(fluidPositions, j)
                {
                    const vector& rJ = fluidPositions[j];
                    
                    scalar rIJ = mag(com - rJ);
                
                    if(rIJ < count*maxR_)
                    {
                        selectPositions[s].append(rJ);
                        selectIds[s].append(j);
                    }
                }
                
                com = vector::zero;
                
                forAll(selectPositions[s], j)
                {
                    com += selectPositions[s][j];
                }
                
                label newMolsCount = selectPositions[s].size();
                
                com /= scalar(newMolsCount);
                
//                     Info << "com = " << com << endl;
            
                count += 1.0;
                
                if(newMolsCount == molsCount)
                {
                    finish = true;
                }
                else
                {
                    selectPositions[s].clear();
                    selectIds[s].clear();                    
                    molsCount = newMolsCount;
                }
            }
            
            // set remaining molecules
            forAll(selectIds[s], i)
            {
                molsChosen[selectIds[s][i]] = s;
            }
        }
        
        
        Info << "sorted lists = " << endl;
        
        forAll(selectPositions, s)
        {
            Info << selectPositions[s].size() << endl;
        }
        
        // reduce lists
        label noOfSegments = 0;
        
        forAll(selectPositions, s)
        {
            if(selectPositions[s].size() > val_)
            {
                noOfSegments++;
            }
        }
        
        List<DynamicList<vector> > sortedPositions(noOfSegments);
        List<DynamicList<label> > sortedIds(noOfSegments);                
            
        label cS = 0;
        
        forAll(selectPositions, s)
        {
            if(selectPositions[s].size() > val_)
            {
                sortedPositions[cS].transfer(selectPositions[s]);
                sortedIds[cS].transfer(selectIds[s]);
                cS++;
            }
        }
        
        
        dropletMolPositions_.clear();
        dropletMolIds_.clear();
        
        if(noOfSegments == 2)
        {
            vector com1 = vector::zero;
            
            forAll(sortedPositions[0], i)
            {
                com1 += sortedPositions[0][i];
            }

            com1 /= sortedPositions[0].size();
            
            vector com2 = vector::zero;
            
            forAll(sortedPositions[1], i)
            {
                com2 += sortedPositions[1][i];
            } 
            
            com2 /= sortedPositions[1].size();
            
            vector c21 = com1 - com2;
            
            scalar sgn = sign(translationVector_ & c21);

            forAll(sortedPositions[0], i)
            {
                dropletMolPositions_.append(sortedPositions[0][i]);
                dropletMolIds_.append(sortedIds[0][i]);
            }
            
            forAll(sortedPositions[1], i)
            {
                dropletMolPositions_.append(sortedPositions[1][i] + sgn*translationVector_);
                dropletMolIds_.append(sortedIds[1][i]);
            }                
        }
        else if(noOfSegments > 2)
        {
            FatalErrorIn("void dropletInfo::calculateField()")
                << "Did not code for more than 2 segements"
                << exit(FatalError);
        }
    }
}

void dropletInfo::calculateField()
{
    
    Info << " Calculate Field " << endl;
    
    if(
        (time_.timeOutputValue() >= startTime_) // skipping some known issues in first time step
        && (time_.timeOutputValue() <= endTime_)
    )
    {
        readFromStorage();
        reconstructDroplet();
//         measureContactAngle();
        
    }
}
    
void dropletInfo::writeField()
{
    // write xmol file
    {
        const reducedUnits& rU = rU_;

        fileName fName(outputPath_/time_.timeName()+"_shifted_"+fieldName_+".xmol");
    
        OFstream str(fName);
    
        str << dropletMolPositions_.size() << nl << "polyMoleculeCloud site positions in angstroms" << nl;

        // for all processors
        forAll(dropletMolPositions_, i)
        {
            vector rI = dropletMolPositions_[i]*rU.refLength()*1e10;

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
