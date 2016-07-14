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
    agentConfiguration

Description


SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef agentConfiguration_H
#define agentConfiguration_H

#include "IOdictionary.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "agent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agentConfiguration Declaration
\*---------------------------------------------------------------------------*/

class agentConfiguration
{

protected:

    // Protected data

        const fvMesh& mesh_;

        agentCloud& cloud_;

        const dictionary& initialiseDict_;

        label nAgentsAdded_;

    // Private Member Functions
    
        // insert knowing from beforehand what cell info is
        void insertAgent
        (
            const point& position,
            const label cell,
            const label tetFace,
            const label tetPt,
            const label& id,          
            const bool& frozen,
            const vector& v
        );
        
//         // insert and determine cell info in function
//         void insertMolecule
//         (
//             const point& position,
//             const label& id,         
//             const bool& tethered,
//             const bool& frozen,
//             const scalar& temperature,
//             const vector& bulkVelocity
//         );
        

        void deleteAgent
        (
            agent& mol
        );        
        
        vector equipartitionLinearVelocity
        (
            scalar temperature,
            scalar mass
        );

public:

    //- Runtime type information
    TypeName("agentConfiguration");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            agentConfiguration,
            dictionary,
            (
                agentCloud& cloud,
                const dictionary& dict
            ),
            (cloud, dict)
        );
   
    // Constructors

        //- Construct from components
        agentConfiguration
        (
            agentCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<agentConfiguration> New
        (
            agentCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~agentConfiguration();


    // Member Functions

        virtual void setInitialConfiguration() = 0;

        const label& nAgents() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
