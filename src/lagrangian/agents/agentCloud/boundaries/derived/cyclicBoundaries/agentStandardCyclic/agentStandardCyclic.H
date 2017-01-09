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
    agentStandardCyclic

Description


SourceFiles
    agentStandardCyclic.C


\*---------------------------------------------------------------------------*/

#ifndef agentStandardCyclic_H
#define agentStandardCyclic_H

#include "agentCyclicBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "agentCloud.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agentStandardCyclic Declaration
\*---------------------------------------------------------------------------*/

class agentStandardCyclic
:
	public agentCyclicBoundary
{

private:

    // Private data
        dictionary propsDict_;



public:

    //- Runtime type information
    TypeName("agentStandardCyclic");

    // Constructors

        //- Construct from components
        agentStandardCyclic
        (
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~agentStandardCyclic();


    // Member Functions

        void calculateProperties();

        //- create an initial configuration
        void initialConfiguration();

        //- control the molecules
        void control
        (
            agent& p,
            const cyclicPolyPatch& cpp,
            agent::trackingData& td
            
//             const label& patchi
        );



        //- output data
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );





        // Access





};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
