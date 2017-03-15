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
    hogmanayControl

Description
    Add agents to a box based on a known flux 
    Not parallelised

SourceFiles
    hogmanayControlI.H
    hogmanayControl.C
    hogmanayControlIO.C

\*---------------------------------------------------------------------------*/

#ifndef hogmanayControl_H
#define hogmanayControl_H

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
                           Class hogmanayControl Declaration
\*---------------------------------------------------------------------------*/

class hogmanayControl
:
	public agentController
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        //- ids of agents
        List<label> topAgentIds_;
        List<label> bottomAgentIds_;

        List<label> rightAgentIds_;
        List<label> leftAgentIds_;                
        
        boundedBox boxBottom_;
        boundedBox boxTop_;
        
        scalar meanDesSpeedE_;
        scalar desSpeedRangeE_;       

        scalar meanDesSpeedF_;
        scalar desSpeedRangeF_;   
        
        // will force 
        
        scalar tau_;

        scalar frac1_;
        scalar frac2_;
        scalar frac3_;            
        
        
        // borders 
        
        List<vectorList> borderList_;
        List<DynamicList<label> > interactionList_;  
        
        List<vectorList> midPoints_;
        List<vectorList> normalVectors_;
        
        vectorList centrePoints_;
        
        autoPtr<agentWallForce> attractiveForceModel_;
//         autoPtr<agentWallForce> repulsiveForceModel_;
        
        // distance at which will force is on/off
//         scalar X1_;
//         scalar X2_;
//         scalar X3_;
        
        scalar Y_;
        scalar YMax_;
        scalar XRight_;
        scalar XLeft_;
        
        // searchTime 
        scalar fireworkTime_;
        scalar evacuateTime_;
//         scalar deltaT_;
        
        // throwing boxes 
        
//         List<boundedBox> boxes_;
        
    // Private Member Functions

//         void readProperties();
        
        void setBoundBox
        (
            const dictionary& propsDict,
            boundedBox& bb,
            const word& name 
        );
     
        
        scalar gaussianDistribution
        (
            const scalar& mean,
            const scalar& range
        );
        
        void initialiseBorders();


        
public:

    //- Runtime type information
    TypeName("hogmanayControl");

    // Constructors

        //- Construct from components
        hogmanayControl
        (
            Time& t,
            agentCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~hogmanayControl();


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
