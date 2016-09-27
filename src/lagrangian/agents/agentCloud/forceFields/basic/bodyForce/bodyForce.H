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
    bodyForce

Description
    Basic/abstract class of a state controller


SourceFiles
    bodyForceI.H
    bodyForce.C
    bodyForceIO.C

\*---------------------------------------------------------------------------*/

#ifndef bodyForce_H
#define bodyForce_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "agent.H"
#include "writeTimeData.H"
#include "selectAgentIds.H"
#include "reducedUnits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class bodyForce Declaration
\*---------------------------------------------------------------------------*/

class bodyForce
{

protected:

    // Protected data

        const fvMesh& mesh_;

        agentCloud& cloud_;

        Time& time_;



    // Private Member Functions

public:

    //- Runtime type information
    TypeName("bodyForce");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            bodyForce,
            dictionary,
            (
                agentCloud& cloud,
                Time& t,
                const dictionary& dict
            ),
            (cloud, t, dict)
        );
   
    // Constructors

        //- Construct from components
        bodyForce
        (
            agentCloud& cloud,
            Time& t,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<bodyForce> New
        (
            agentCloud& cloud,
            Time& t,
            const dictionary& dict
        );


    // Destructor

        virtual ~bodyForce();


    // Member Functions


        //- create an initial configuration
        virtual void initialConfiguration() = 0;
        
        // inside loop of all agents
        virtual void force(agent* p) = 0;
        
        // new force definition; requires loop over all agents
        virtual void newForce() = 0;                
        
        
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
