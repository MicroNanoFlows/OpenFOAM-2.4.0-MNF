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
    agentCyclicBoundary

Description

    Used for applying instantaneous actions on molecules colliding with a cyclic 
    boundary.

SourceFiles
    agentCyclicBoundaryI.H
    agentCyclicBoundary.C
    agentCyclicBoundaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef agentCyclicBoundary_H
#define agentCyclicBoundary_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "polyMesh.H"
#include "agent.H"
#include "Time.H"
#include "writeTimeData.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agentCyclicBoundary Declaration
\*---------------------------------------------------------------------------*/

class agentCyclicBoundary
{

protected:

    // Protected data

        const fvMesh& mesh_;

        agentCloud& cloud_;

        Time& time_;

        //- name of this patch
        word patchName_;
        
        //- id on boundaryMesh
        label patchId_;
        
        //- writing out properties
        bool writeInTimeDir_;
        bool writeInCase_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("agentCyclicBoundary");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            agentCyclicBoundary,
            dictionary,
            (
                Time& t,
                const polyMesh& mesh,
                agentCloud& cloud,
                const dictionary& dict
            ),
            (t, mesh, cloud, dict)
        );
   
    // Constructors

        //- Construct from components
        agentCyclicBoundary
        (
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<agentCyclicBoundary> New
        (
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~agentCyclicBoundary();


    // Member Functions


        virtual void calculateProperties() = 0;

        //- create an initial configuration
        virtual void initialConfiguration() = 0;

        //- control the molecule
        virtual void control
        (
            agent& p,
            const cyclicPolyPatch& cpp,
            agent::trackingData& td
        ) = 0;


        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;
        
        // Access

            const word& patchName() const;
            const label& patchId() const;
            const bool& writeInTimeDir() const;
            const bool& writeInCase() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
