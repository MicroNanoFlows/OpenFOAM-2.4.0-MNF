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
    agentWillForceScheduled

Description


SourceFiles
    agentWillForceScheduledI.H
    agentWillForceScheduled.C
    agentWillForceScheduledIO.C

\*---------------------------------------------------------------------------*/

#ifndef agentWillForceScheduled_H
#define agentWillForceScheduled_H

#include "bodyForce.H"
#include "vector.H"
#include "polyMesh.H"
#include "agentCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agentWillForceScheduled Declaration
\*---------------------------------------------------------------------------*/

class agentWillForceScheduled
:
	public bodyForce
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        //- ids of agents
        List<label> agentIds_;
        
//         scalar desiredSpeed_;
        
        vector desiredDirection_;
        
        scalar tau_;

//         scalar stdev_;
        
        // size = number of agents
//         List<scalar> agentVelocities_;
        
    // Private Member Functions


public:

    //- Runtime type information
    TypeName("agentWillForceScheduled");

    // Constructors

        //- Construct from components
        agentWillForceScheduled
        (
            agentCloud& cloud,
            Time& t,
            const dictionary& dict
        );

    // Destructor

        ~agentWillForceScheduled();


    // Member Functions
        
        //- create an initial configuration
        void initialConfiguration();
        
        // inside loop of all agents
        void force(agent* p);
        
        // new force definition; requires loop over all agents
        void newForce();

        // Access
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
