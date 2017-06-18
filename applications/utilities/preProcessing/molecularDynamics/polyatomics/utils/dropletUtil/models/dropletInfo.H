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
    dropletInfo

Description

    
    
    
SourceFiles
    dropletInfo.C


\*---------------------------------------------------------------------------*/

#ifndef dropletInfo_H
#define dropletInfo_H

#include "utilField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "distribution.H"
#include "writeTimeData.H"
#include "boundedBox.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dropletInfo Declaration
\*---------------------------------------------------------------------------*/

class dropletInfo
:
    public utilField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        word fieldName_;

        fileName nameFile1_;
        fileName nameFile2_;        
        fileName nameFile3_;           
        fileName nameFile4_; 
        fileName nameFile5_; 
 
        word translate_;

        vectorField positions_;
        scalarField energies_;
        vectorField forces_;
        scalarField nPairs_;
        List<bool> fluidMols_; 

        List<bool> newMolsToPick_;
        
        label nFluidMols_;
        
        label nSelectedMols_;
        
        scalar startTime_;
        scalar endTime_;
        
       
        scalar deltaT_;

        scalar maxR_;
        
        label nSegs_;

        label val_;
        
        vector translationVector_;
        
        DynamicList<vector> dropletMolPositions_;
        DynamicList<label> dropletMolIds_;
        
    // Private Member Functions
        
        void readFromStorage();
        void reconstructDroplet();
            
        

 public:

    //- Runtime type information
    TypeName("dropletInfo");

    // Constructors

        //- Construct from components
        dropletInfo
        (
            Time& t,
            const polyMesh& mesh,
            const reducedUnits& rU,
            const dictionary& dict
        );

    // Destructor

        ~dropletInfo();


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
