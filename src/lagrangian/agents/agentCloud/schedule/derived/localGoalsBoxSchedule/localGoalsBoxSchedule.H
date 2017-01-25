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
    localGoalsBoxSchedule

Description


SourceFiles
    localGoalsBoxScheduleI.H
    localGoalsBoxSchedule.C
    localGoalsBoxScheduleIO.C

\*---------------------------------------------------------------------------*/

#ifndef localGoalsBoxSchedule_H
#define localGoalsBoxSchedule_H

#include "scheduleModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "agentCloud.H"
#include "boundedBox.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class localGoalsBoxSchedule Declaration
\*---------------------------------------------------------------------------*/

class localGoalsBoxSchedule
:
	public scheduleModel
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        //- ids of agents
        List<label> agentIds_;
        
        List<boundedBox> boxes_;
        // list of centre points of the box
        vectorList destinations_;
//         scalarList regionRadii_;
        label nScheduledDestinations_;
        
        scalar timeIndex_;
        label counter_;
        scalar deltaT_;
        
    // Private Member Functions
        
        void setBoundBoxes();


public:

    //- Runtime type information
    TypeName("localGoalsBoxSchedule");

    // Constructors

        //- Construct from components
        localGoalsBoxSchedule
        (
            Time& time,
            agentCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~localGoalsBoxSchedule();

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
