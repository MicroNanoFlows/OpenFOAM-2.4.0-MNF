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
    localGoalsSchedule

Description


SourceFiles
    localGoalsScheduleI.H
    localGoalsSchedule.C
    localGoalsScheduleIO.C

\*---------------------------------------------------------------------------*/

#ifndef localGoalsSchedule_H
#define localGoalsSchedule_H

#include "scheduleModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "agentCloud.H"
#include "boundedBox.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class localGoalsSchedule Declaration
\*---------------------------------------------------------------------------*/

class localGoalsSchedule
:
	public scheduleModel
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        //- ids of agents
        List<label> agentIds_;
        
        vectorList destinations_;
        scalarList regionRadii_;
        scalarList timeAtDestination_;
        
        label nScheduledDestinations_;
        
        scalar timeIndex_;
        label counter_;
        scalar deltaT_;
        
    // Private Member Functions


public:

    //- Runtime type information
    TypeName("localGoalsSchedule");

    // Constructors

        //- Construct from components
        localGoalsSchedule
        (
            Time& time,
            agentCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~localGoalsSchedule();

    // Member Functions
        
        void initialConfiguration();
        
        // is called every time-step to set new destinations of agents
        void setSchedule();

                

        // Access
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
