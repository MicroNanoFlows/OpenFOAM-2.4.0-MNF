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
    agentPatchBoundary

Description

    Used for applying instantaneous actions on molecules colliding with a boundary
    patch.

SourceFiles
    agentPatchBoundaryI.H
    agentPatchBoundary.C
    agentPatchBoundaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef agentPatchBoundary_H
#define agentPatchBoundary_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "agent.H"
#include "writeTimeData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agentPatchBoundary Declaration
\*---------------------------------------------------------------------------*/

class agentPatchBoundary
{

protected:

    // Protected data

        const fvMesh& mesh_;

        Time& t_;

        agentCloud& cloud_;

        //- subDictionary containing the properties
//         dictionary boundaryDict_;

        //- name of boundary patch
        word patchName_;
        label patchId_;

        //- faces on the boundary
        labelList faces_;
        label nFaces_;

        //- cells corresponding to the boundary faces
        labelList cells_;




            //- writing out properties
            bool writeInTimeDir_;
            bool writeInCase_;


    // Private Member Functions


public:

    //- Runtime type information
    TypeName("agentPatchBoundary");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            agentPatchBoundary,
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
        agentPatchBoundary
        (
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<agentPatchBoundary> New
        (
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~agentPatchBoundary();


    // Member Functions

        void setBoundaryCells();

        //- calculate any required properties
        virtual void calculateProperties() = 0;

        virtual void initialConfiguration() = 0;

        //- control the agent
        virtual void controlMol
        (
            agent& mol,
            agent::trackingData& td
        ) = 0;

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        //- update properties from a modified dictionary
//         virtual void updateProperties(const dictionary&) = 0;

        void updateBoundaryProperties(const dictionary&);

        // Access

            const labelList& controlPatch() const;

            //- return the adjacent cells
            const labelList& controlZone() const;

            //- return the name of the control patch
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
