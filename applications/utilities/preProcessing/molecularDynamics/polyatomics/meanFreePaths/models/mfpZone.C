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

#include "mfpZone.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(mfpZone, 0);

addToRunTimeSelectionTable(mfpField, mfpZone, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
mfpZone::mfpZone
(
    Time& t,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mfpField(t, mesh, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.lookup("fieldName")),
    nameFile1_(propsDict_.lookup("nameFile_collectorOfFreePaths")),
    nameFile2_(propsDict_.lookup("nameFile_freePaths")),
    nameFile3_(propsDict_.lookup("nameFile_nCollisions")),
    nameFile4_(propsDict_.lookup("nameFile_startCollisionPositions_X")),
    nameFile5_(propsDict_.lookup("nameFile_startCollisionPositions_Y")),
    nameFile6_(propsDict_.lookup("nameFile_startCollisionPositions_Z")),
    nameFile7_(propsDict_.lookup("nameFile_endCollisionPositions_X")),
    nameFile8_(propsDict_.lookup("nameFile_endCollisionPositions_Y")),
    nameFile9_(propsDict_.lookup("nameFile_endCollisionPositions_Z")),
    nameFile10_(propsDict_.lookup("nameFile_startCollisionTimes")),
    nameFile11_(propsDict_.lookup("nameFile_endCollisionTimes"))
{
    

    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mfpZone::~mfpZone()
{}





// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mfpZone::createField()
{
    readFromStorage();
    
    //test 
    forAll(nCollisions_, i)
    {
        if(nCollisions_[i] < 1)
        {
            FatalErrorIn("void mfpZone::createField()")
                << "File: " << nameFile3_
                << ", contains an entry which has no collisions" 
                << nl << "   -> this will result in a bug in this code, "
                << " so I've decided to terminate until solved."
                << abort(FatalError);            
        }
    }
    
    nMols_ = nCollisions_.size();    
}

// void mfpZone::setBoundBox
// (
//     const dictionary& propsDict,
//     boundedBox& bb,
//     const word& name 
// )
// {
//     const dictionary& dict(propsDict.subDict(name));
//     
//     vector startPoint = dict.lookup("startPoint");
//     vector endPoint = dict.lookup("endPoint");
// 
//     bb.resetBoundedBox(startPoint, endPoint);
// }

void mfpZone::calculateField()
{
    // find the last free path that made a fresh start 
    label index = -1;
    scalar tStart = 0.0;
    
    for (label i = 0; i < nMols_; i++)
    {
        if(startCollTimes_[i][0] > tStart)
        {
            tStart = startCollTimes_[i][0];
            index = i;
        }
    }
    
    Info << "start at time = " << tStart << ", mol = " << index << endl;
}
    
void mfpZone::writeField()
{

//     scalarField timeField (mfpZone_.size(), 0.0);
//     
//     mfpZone.transfer(mfpZone_);
//     mfpZone_.clear();
// 
//     const scalar& deltaT = time_.time().deltaT().value();
//     
//     forAll(timeField, i)
//     {
//         timeField[timeField.size()-i-1]=runTime.timeOutputValue()-(deltaT*i);
//     }            
//     
//     writeTimeData
//     (
//         outputPath_,
//         "mfpZone_"+fieldName_+".xyz",
//         timeField,
//         mfpZone,
//         true
//     );

        
       
}

void mfpZone::readFromStorage()
{
    Info << "reading from storage" << endl;
    
    fileName pathName(inputPath_);
       
    {
        IFstream file(pathName/nameFile1_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> collectorOfFreePaths_;
        }
    }
    {
        IFstream file(pathName/nameFile2_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> freePaths_;
        }
    }    
    {
        IFstream file(pathName/nameFile3_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> nCollisions_;
        }
    }
    {
        IFstream file(pathName/nameFile4_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> startCollPosX_;
        }
    }
    {
        IFstream file(pathName/nameFile5_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> startCollPosY_;
        }
    }
    {
        IFstream file(pathName/nameFile6_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> startCollPosZ_;
        }
    }
    {
        IFstream file(pathName/nameFile7_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> endCollPosX_;
        }
    }
    {
        IFstream file(pathName/nameFile8_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> endCollPosY_;
        }
    }
    {
        IFstream file(pathName/nameFile9_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> endCollPosZ_;
        }
    }    
    {
        IFstream file(pathName/nameFile10_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> startCollTimes_;
        }
    }
    {
        IFstream file(pathName/nameFile11_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> endCollTimes_;
        }
    } 
    
//     Info << "collectorOfFreePaths = " << collectorOfFreePaths_ << endl;
//     Info << "freePaths = " << freePaths_ << endl;
//     Info << "nCollisions = " << nCollisions_ << endl;        
          
    
}

} // End namespace Foam

// ************************************************************************* //
