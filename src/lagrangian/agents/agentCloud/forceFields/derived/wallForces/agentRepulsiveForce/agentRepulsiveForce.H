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
    agentRepulsiveForceI.H
    agentRepulsiveForce.C
    agentRepulsiveForceIO.C

\*---------------------------------------------------------------------------*/

#ifndef agentRepulsiveForce_H
#define agentRepulsiveForce_H


#include "agentWallForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agentRepulsiveForce Declaration
\*---------------------------------------------------------------------------*/

class agentRepulsiveForce
:
    public agentWallForce
{

private:

    // Private data

        dictionary propsDict_;

        vector force_;


public:

    //- Runtime type information
    TypeName("agentRepulsiveForce");

    // Static data members

    // Constructors

        //- construct from components
        agentRepulsiveForce
        (
            Time& time,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~agentRepulsiveForce();

    // Member Functions

//         void initialConfig
//         (
//             const scalar& density,
//             const radialDistribution& rdf,
//             const moleculeCloud& mD
//         );

        //- output data
//         void output(const Time& runTime, const polyMesh& mesh);


        // Access

        // Check

        // Edit

            void updateForce();

            vector force(const vector& position);
            vector force(const scalar& time);
            vector force(const scalar& time, const vector& position);
            
        // Write


            void write
            (
                const fileName& fixedPathName,
                const fileName& timePath
            );


//             void updateProperties
//             (
//                 const dictionary& dict
//             );

    // Member Operators

//         void operator=(const agentRepulsiveForce&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const agentRepulsiveForce&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
