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

Description

\*---------------------------------------------------------------------------*/

#include "laneFormationControl.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(laneFormationControl, 0);

addToRunTimeSelectionTable(agentController, laneFormationControl, dictionary);




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
laneFormationControl::laneFormationControl
(
    Time& t,
    agentCloud& cloud,
    const dictionary& dict
)
:
    agentController(t,  cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
//     agentId_()
{
    writeInTimeDir_ = false;
    writeInCase_ = false;

    {
        leftAgentIds_.clear();

        selectAgentIds ids
        (
            cloud_.cP(),
            propsDict_,
            "leftAgentIds"
        );

        leftAgentIds_ = ids.agentIds();
    }    
    
    {
        rightAgentIds_.clear();

        selectAgentIds ids
        (
            cloud_.cP(),
            propsDict_,
            "rightAgentIds"
        );

        rightAgentIds_ = ids.agentIds();
    }
    
    tau_ = readScalar(propsDict_.lookup("tau"));
    
    frac1_ = readScalar(propsDict_.lookup("frac1"));

        
    // borders 
    leftDir_ = propsDict_.lookup("leftDir");
    rightDir_ = propsDict_.lookup("rightDir"); 
    
    
    leftDir_ /= mag(leftDir_);
    rightDir_ /= mag(rightDir_);
    
//     X1_ = readScalar(propsDict_.lookup("X1"));    
//     X2_ = readScalar(propsDict_.lookup("X2")); 
//     X3_ = readScalar(propsDict_.lookup("X3")); 
//     
    
//     Y_ = readScalar(propsDict_.lookup("Y")); 
    
 
    
//     deltaT_ = time_.deltaT().value();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

laneFormationControl::~laneFormationControl()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void laneFormationControl::initialConfiguration()
{}

void laneFormationControl::controlBeforeVelocityI()
{}

void laneFormationControl::controlBeforeMove()
{}

void laneFormationControl::controlBeforeForces()
{

}

void laneFormationControl::controlDuringForces
(
    agent* molI,
    agent* molJ
)
{}

void laneFormationControl::controlAfterForces()
{
    {
        IDLList<agent>::iterator p(cloud_.begin());
        
        for (p = cloud_.begin(); p != cloud_.end(); ++p)
        {
            if(p().eventTracker() == 0)
            {
                if(findIndex(leftAgentIds_, p().id()) != -1)
                {
                    vector desDir = leftDir_;
                    vector force = (p().desiredSpeed()*desDir - p().v())*p().mass() / tau_; 
                    p().f() += force;
                    vector R = vector(0.0, 2.0*cloud_.rndGen().scalar01()-1.0, 0.0);
                    R /= mag(R);
                    p().f() += R*frac1_*mag(force);       

                }
                else if (findIndex(rightAgentIds_, p().id()) != -1)
                {
                    vector desDir = rightDir_;
                    vector force = (p().desiredSpeed()*desDir - p().v())*p().mass() / tau_; 
                    p().f() += force;                    
                    vector R = vector(0.0, 2.0*cloud_.rndGen().scalar01()-1.0, 0.0);
                    R /= mag(R);
                    p().f() += R*frac1_*mag(force);       
                }
            }
        }
                
    }
}




void laneFormationControl::controlAfterVelocityII()
{}

void laneFormationControl::calculateProperties()
{}

void laneFormationControl::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
    }
}


} // End namespace Foam

// ************************************************************************* //
