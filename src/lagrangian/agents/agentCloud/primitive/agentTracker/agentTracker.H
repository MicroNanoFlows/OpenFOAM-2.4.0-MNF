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
    Foam::agentTracker

Description



SourceFiles
    agentTrackerI.H
    agentTracker.C
    agentTrackerIO.C

\*---------------------------------------------------------------------------*/

#ifndef agentTracker_H
#define agentTracker_H

#include "scalar.H"
#include "IOdictionary.H"
#include "polyMesh.H"
#include "Time.H"
#include "agent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class agentTracker Declaration
\*---------------------------------------------------------------------------*/

class agentTracker
{
    // Private data

        agentCloud& cloud_;
    
        DynamicList<agent*> agentsOnProc_;
        
        DynamicList<label> agentTrackingNumbers_;
        

        // tool variables
        
//         label trackingNumber_;
        
        label id_;
        

    // Private Member Functions



public:
    
        friend class agent;


    // Constructors


        //- Construct from components
        agentTracker
        (
            agentCloud& cloud
        );


    // Destructor

        ~agentTracker();


    // Member Functions

            void track();
        
        // Access
            
            // this is a two step function
            
            // call this first (to avoid errors)
            inline bool isTrackingNumberAvailable
            (
                const label& trackingNumber
            );

            // then call this function to get the agent
            inline agent* getAgent();            
            


    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const agentTracker&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "agentTrackerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
