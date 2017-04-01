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
    meanFreePath

Description
    This tool works well if you initialise a set of molecules with tracking Numbers
    that range between 0 and N. In the case of doing complex things,  like using
    FADE, I would recommend you reinitiliase the trackingNumbers.

SourceFiles
    meanFreePathI.H
    meanFreePath.C
    meanFreePathIO.C

\*---------------------------------------------------------------------------*/

#ifndef meanFreePath_H
#define meanFreePath_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "distribution.H"
#include "boundedBox.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class meanFreePath Declaration
\*---------------------------------------------------------------------------*/

class meanFreePath
:
	public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        propertyField fields_;

        word fieldName_;

        List<label> molIds_;
        
        fileName nameFile1_;
        fileName nameFile2_;        
        fileName nameFile3_;           
        
        // collision distance
        scalar dCol_;
        
        // MD timestep
        scalar deltaT_;
        
        label nMols_;
        
        label listSize_;
        
        List<DynamicList<scalar> > collectorOfFreePaths_;        
        List<scalar> freePaths_;
        List<scalar> nCollisions_;        
        List<scalar> Rold_;

        bool testInitialTimeScale_;
        
    // Private Member Functions

        void setRolds();
        
        void readFromStorage();
        
        void writeToStorage();
        
public:

    //- Runtime type information
    TypeName("meanFreePath");

    // Constructors

        //- Construct from components
        meanFreePath
        (
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );

    // Destructor

        ~meanFreePath();


    // Member Functions

        //- create initial field
        void createField();

        void afterForce();
        
        //- calculate field
        void calculateField();

        void measureDuringForceComputation
        (
            polyMolecule* molI,
            polyMolecule* molJ
        );

        void measureDuringForceComputationSite
        (
            polyMolecule* molI,
            polyMolecule* molJ,
            label sI,
            label sJ
        );

        const propertyField& fields() const;

        //- write field
        void writeField();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
