/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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


Description

SourceFiles
    uniform3DBins.C

\*---------------------------------------------------------------------------*/

#ifndef uniform3DBins_H
#define uniform3DBins_H


#include "threeDimBinModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uniform3DBins Declaration
\*---------------------------------------------------------------------------*/

class uniform3DBins
:
    public threeDimBinModel
{

private:

    // Private data

        dictionary propsDict_;

        vector startPoint_;
        vector startPointX_;
        vector startPointY_;
        vector startPointZ_;
        vector endPoint_;
        vector endPointX_;
        vector endPointY_;
        vector endPointZ_;
        scalar rSEMagX_;
        scalar rSEMagY_;
        scalar rSEMagZ_;


        vector nBins_;
        vector binWidth_;
        scalar area_;

public:

    //- Runtime type information
    TypeName("uniform3DBins");

    // Static data members

    // Constructors

        //- construct from components
        uniform3DBins
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

    // Selectors

    // Destructor

        ~uniform3DBins();

    // Member Functions

        // Access

        // Check

        // Edit

            List<label> isPointWithinBin
            (
                const vector& rI,
                const label& cellI
            );

            vectorField binPositionsX();
            scalarField binPositionsY();
            scalarField binPositionsZ();

            vectorField binPositionsXYZ();

            vector position(const vector& h, const scalar& r, const scalar& theta);

            // 0 length, 1 radius, 2 angle
            List<label> nBins();

            scalar binVolume(const label& n);

            void write
            (
                const fileName& path,
                const word& name
            );

            

        // Write

    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
