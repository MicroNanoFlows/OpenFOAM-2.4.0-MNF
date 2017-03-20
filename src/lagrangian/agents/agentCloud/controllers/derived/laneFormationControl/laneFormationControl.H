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
    laneFormationControl

Description
    Add agents to a box based on a known flux 
    Not parallelised

SourceFiles
    laneFormationControlI.H
    laneFormationControl.C
    laneFormationControlIO.C

\*---------------------------------------------------------------------------*/

#ifndef laneFormationControl_H
#define laneFormationControl_H

#include "agentController.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "agentCloud.H"
#include "boundedBox.H"
#include "agentWallForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class laneFormationControl Declaration
\*---------------------------------------------------------------------------*/

class laneFormationControl
:
	public agentController
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        //- ids of agents
        List<label> leftAgentIds_;
        List<label> rightAgentIds_;


        
        // will force 
        
        scalar tau_;

        scalar frac1_;
         
        vector leftDir_;
        vector rightDir_;
        
        
    // Private Member Functions


        
public:

    //- Runtime type information
    TypeName("laneFormationControl");

    // Constructors

        //- Construct from components
        laneFormationControl
        (
            Time& t,
            agentCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~laneFormationControl();


    // Member Functions
    // Member Functions
        
        //- create an initial configuration
        void initialConfiguration();

        //- control molecules at various stages of time-step
        void controlBeforeVelocityI();
        void controlBeforeMove();
        void controlBeforeForces();
        void controlDuringForces
        (
            agent* molI,
            agent* molJ
        );        
        void controlAfterForces();
        void controlAfterVelocityII();
        
        //- calculate any required properties (end of time-step)
        void calculateProperties(); 
        
        //- output data (optional)
        void output
        (
            const fileName& fixedPathName, 
            const fileName& timePath
        );
        
        //- update properties from a modified dictionary (optional)
//         void updateProperties(const dictionary&);

        // Access
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
