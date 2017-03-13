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
    jamaratControl

Description
    Add agents to a box based on a known flux 
    Not parallelised

SourceFiles
    jamaratControlI.H
    jamaratControl.C
    jamaratControlIO.C

\*---------------------------------------------------------------------------*/

#ifndef jamaratControl_H
#define jamaratControl_H

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
                           Class jamaratControl Declaration
\*---------------------------------------------------------------------------*/

class jamaratControl
:
	public agentController
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        //- ids of agents
        List<label> agentIds_;

        
        //- id of inserted agent
        label agentId_;

        boundedBox samplingBox_;
        boundedBox controlBox_;
        
        label N_;
        
        scalar meanMass_;
        scalar massRange_;
        scalar meanRadius_;
        scalar radiusRange_;
        scalar meanDesSpeed_;
        scalar desSpeedRange_;
        
        
        // will force 
        
        vector desiredDirection_;
        scalar tau_;
        scalar fracP_;
        scalar fracR_;
        // borders 
        
        List<vectorList> borderList_;
        List<DynamicList<label> > interactionList_;  
        
        List<vectorList> midPoints_;
        List<vectorList> normalVectors_;
        
        vectorList centrePoints_;
        
        List<vector> sideTop_;
        List<vector> sideBottom_;
        
        
        autoPtr<agentWallForce> attractiveForceModel_;
        autoPtr<agentWallForce> repulsiveForceModel_;
        
        // distance at which will force is on/off
        scalar X1_;
        scalar X2_;
        scalar X3_;
        
        // searchTime 
        scalar searchTime_;
        scalar throwingTime_;
        scalar deltaT_;
        
        // throwing boxes 
        
        List<boundedBox> boxes_;
        
    // Private Member Functions

//         void readProperties();
        
        void setBoundBox
        (
            const dictionary& propsDict,
            boundedBox& bb,
            const word& name 
        );
     
        void insertAgents(const label& N);
        
        scalar gaussianDistribution
        (
            const scalar& mean,
            const scalar& range
        );
        
        void initialiseBorders();
        
        void setBoundBoxes();
        
//         bool isPointWithinBorder(label index, const vector& r);
        
        void reflect(label index, agent* p);
        
public:

    //- Runtime type information
    TypeName("jamaratControl");

    // Constructors

        //- Construct from components
        jamaratControl
        (
            Time& t,
            agentCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~jamaratControl();


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
