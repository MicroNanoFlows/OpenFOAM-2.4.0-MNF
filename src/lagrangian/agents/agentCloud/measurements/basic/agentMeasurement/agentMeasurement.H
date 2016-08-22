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
    agentMeasurement

Description

SourceFiles
    agentMeasurementI.H
    agentMeasurement.C
    agentMeasurementIO.C

\*---------------------------------------------------------------------------*/

#ifndef agentMeasurement_H
#define agentMeasurement_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "agent.H"
#include "writeTimeData.H"
#include "selectAgentIds.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agentMeasurement Declaration
\*---------------------------------------------------------------------------*/

class agentMeasurement
{

protected:

    // Protected data

        const fvMesh& mesh_;

        agentCloud& cloud_;

        Time& time_;

        fileName casePath_;

        fileName timePath_;

        bool measureInterForces_;
        bool measureInterForcesSites_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("agentMeasurement");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            agentMeasurement,
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
        agentMeasurement
        (
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<agentMeasurement> New
        (
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud,
            const dictionary& dict
        );



    // Destructor

        virtual ~agentMeasurement();


    // Member Functions


        const fileName& casePath() const;
        fileName& casePath();

        const fileName& timePath() const;
        fileName& timePath();

        const bool& measureInterForces() const;
        bool& measureInterForces();

        const bool& measureInterForcesSites() const;
        bool& measureInterForcesSites();

        //- create initial field
        virtual void createField() = 0;

        //- calculate field
        virtual void calculateField() = 0;

        //- write field
        virtual void writeField() = 0;

        virtual void measureDuringForceComputation
        (
            agent* molI,
            agent* molJ
        ) = 0;

        virtual void measureDuringForceComputationSite
        (
            agent* molI,
            agent* molJ,
            label sI,
            label sJ
        ) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
