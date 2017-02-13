/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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
    Foam::tiltedCylinder

Description

SourceFiles
    tiltedCylinderI.H
    tiltedCylinder.C
    tiltedCylinderIO.C

\*---------------------------------------------------------------------------*/

#ifndef tiltedCylinder_H
#define tiltedCylinder_H

#include "scalar.H"
#include "IOdictionary.H"
#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class tiltedCylinder Declaration
\*---------------------------------------------------------------------------*/

class tiltedCylinder
{
    // Private data

        vector startPoint_;
        
        vector endPoint_;
        
        scalar radius_;
        
        vector normal_;
        
        scalar magRSE_;

    // Private Member Functions

//         void calcRefValues();

        //- Disallow default bitwise copy construct
//         tiltedCylinder(const tiltedCylinder&);

        //- Disallow default bitwise assignment
//         void operator=(const tiltedCylinder&);


public:

    // Static data members

        //- Static data 
//         static const scalar kb;
//         static const scalar elementaryCharge;
//         static const scalar vacuumPermittivity;

    // Constructors

        //- Construct with no argument
        tiltedCylinder();

        //- Construct from components
        tiltedCylinder
        (
            const vector& startPoint,
            const vector& endPoint,
            const scalar& radius
        );


    // Destructor

        ~tiltedCylinder();


    // Member Functions

            void setBox();
            
            void setBox
            (
                const vector& startPoint,
                const vector& endPoint,
                const scalar& radius                
            );
            

        // Access

            inline bool contains(const vector& p);

 


    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const tiltedCylinder&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tiltedCylinderI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
