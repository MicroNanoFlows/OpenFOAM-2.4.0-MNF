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
    siteI.H
    site.C
    siteIO.C

\*---------------------------------------------------------------------------*/

#ifndef site_H
#define site_H

#include "IFstream.H"
#include "OFstream.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class site Declaration
\*---------------------------------------------------------------------------*/

class site
{
    // Private data

        vector sitePosition_;
        vector siteForce_;

    // Private Member Functions

public:

    // Static data members

    // Constructors

        //- Construct null
        site();

        //- Construct from xs and ys
        site
        (
            const vector& sitePosition,
            const vector& siteForce
        );


    // Selectors

    // Destructor

        ~site();

    // Member Functions

        // Access

            inline const vector& sitePosition() const;
            inline vector& sitePosition();

            inline const vector& siteForce() const;
            inline vector& siteForce();


        // Check

        // Edit

        // Write


    // Member Operators

//         void operator=(const site&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const site&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "siteI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
