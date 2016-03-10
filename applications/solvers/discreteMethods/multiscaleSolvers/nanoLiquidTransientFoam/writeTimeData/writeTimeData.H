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
    writeTimeData

Description

	Stores all the information for the writeTimeData models (state and fluxes),
    defined within the writeTimeDataDict, and selects & builds the models automatically.

\*---------------------------------------------------------------------------*/

#ifndef writeTimeData_H
#define writeTimeData_H

#include "vector.H"
#include "List.H"
#include "IOdictionary.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "volFields.H"
#include "Time.H"
#include "OFstream.H"
#include "IFstream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class writeTimeData Declaration
\*---------------------------------------------------------------------------*/

class writeTimeData
{

    // Private data



public:

    // Constructors

        writeTimeData
        (
            const fileName& pathName,
            const word& fileName,
            const List< Pair<scalar> >& data
        );

        writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const List< Pair<scalar> >& data,
            const label& dummy
        );

        //- scalar field, scalar field
        writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const scalarField& yData
        );

        //- vector field, scalar field, 
        writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const vectorField& xData,
            const scalarField& yData
        );

        //- scalar field, vector field
        writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const vectorField& yData
        );

        //- vector field, vector field
        writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const vectorField& xData,
            const vectorField& yData
        );

        //- scalar field, tensor field
        writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const tensorField& yData
        );

        //- vector field, tensor field
        writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const vectorField& xData,
            const tensorField& yData
        );

        // two scalar fields (with append possible)
        writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const scalarField& yData,
            const bool& dummy
        );

        // one scalar field (with append possible)
        writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const bool& dummy
        );


        // one scalar field one VECTOR field (with append possible)
        writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const vectorField& yData,
            const bool& dummy
        );

        // one VECTOR field (with append possible)
        writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const vectorField& yData,
            const bool& dummy
        );

        // one scalar field one TENSOR field (with append possible)
        writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const tensorField& yData,
            const bool& dummy
        );

        // one TENSOR field (with append possible)
        writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const tensorField& xData,
            const bool& dummy
        );

        ~writeTimeData();


    // Member Functions


    // Ostream Operator

//         friend Ostream& operator<<(Ostream&, const writeTimeData&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
