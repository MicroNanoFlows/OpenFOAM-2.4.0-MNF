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
    mfpField

Description

SourceFiles
    mfpFieldI.H
    mfpField.C
    mfpFieldIO.C

\*---------------------------------------------------------------------------*/

#ifndef mfpField_H
#define mfpField_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "writeTimeData.H"
#include "reducedUnits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mfpField Declaration
\*---------------------------------------------------------------------------*/

class mfpField
{

protected:

    // Protected data

        const fvMesh& mesh_;

        Time& time_;
        
        const reducedUnits& rU_;
        
        fileName inputPath_;

        fileName outputPath_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("mfpField");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            mfpField,
            dictionary,
            (
                Time& t,
                const polyMesh& mesh,
                const reducedUnits& rU,
                const dictionary& dict
            ),
            (t, mesh, rU, dict)
        );
   
    // Constructors

        //- Construct from components
        mfpField
        (
            Time& t,
            const polyMesh& mesh,
            const reducedUnits& rU,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<mfpField> New
        (
            Time& t,
            const polyMesh& mesh,
            const reducedUnits& rU,
            const dictionary& dict
        );



    // Destructor

        virtual ~mfpField();


    // Member Functions


        const fileName& outputPath() const;
        fileName& outputPath();

        const fileName& inputPath() const;
        fileName& inputPath();

        //- create initial field
        virtual void createField() = 0;

        //- calculate field
        virtual void calculateField() = 0;
        
        //- write field
        virtual void writeField() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
