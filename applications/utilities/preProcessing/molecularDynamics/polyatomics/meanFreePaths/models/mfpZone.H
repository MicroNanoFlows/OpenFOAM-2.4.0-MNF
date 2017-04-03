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

Class
    mfpZone

Description

    
    
    
SourceFiles
    mfpZone.C


\*---------------------------------------------------------------------------*/

#ifndef mfpZone_H
#define mfpZone_H

#include "mfpField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
// #include "boundedBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mfpZone Declaration
\*---------------------------------------------------------------------------*/

class mfpZone
:
    public mfpField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

//         boundedBox bb_;
        
        word fieldName_;

        fileName nameFile1_;
        fileName nameFile2_;        
        fileName nameFile3_;           
        fileName nameFile4_; 
        fileName nameFile5_; 
        fileName nameFile6_; 
        fileName nameFile7_; 
        fileName nameFile8_; 
        fileName nameFile9_; 
        fileName nameFile10_; 
        fileName nameFile11_;  

        List<DynamicList<scalar> > collectorOfFreePaths_;
        List<scalar> freePaths_;
        List<scalar> nCollisions_;

        // points at which collisions occur
        List<DynamicList<scalar> > startCollPosX_;
        List<DynamicList<scalar> > startCollPosY_;
        List<DynamicList<scalar> > startCollPosZ_;

        List<DynamicList<scalar> > endCollPosX_;
        List<DynamicList<scalar> > endCollPosY_;
        List<DynamicList<scalar> > endCollPosZ_;

        // times at which collisions occur
        List<DynamicList<scalar> > startCollTimes_;
        List<DynamicList<scalar> > endCollTimes_;
        

        label nMols_;
        
        DynamicList<scalar> mfp_;

        
        
        
    // Private Member Functions
        
        void readFromStorage();

 public:

    //- Runtime type information
    TypeName("mfpZone");

    // Constructors

        //- Construct from components
        mfpZone
        (
            Time& t,
            const polyMesh& mesh,
            const dictionary& dict
        );

    // Destructor

        ~mfpZone();


    // Member Functions

        //- create initial field
        void createField();
        
        //- calculate field
        void calculateField();

        //- write field
        void writeField();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
