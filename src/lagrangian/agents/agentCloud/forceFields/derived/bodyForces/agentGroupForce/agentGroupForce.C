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

#include "agentGroupForce.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agentGroupForce, 0);

addToRunTimeSelectionTable(bodyForce, agentGroupForce, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
agentGroupForce::agentGroupForce
(
    agentCloud& cloud,
    Time& t,
    const dictionary& dict
)
:
    bodyForce(cloud, t, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
//     agentIds_(),
    convergentSpeed_(readScalar(propsDict_.lookup("convergentSpeed"))),
//     desiredDirection_(propsDict_.lookup("desiredDirection")),
    radius_(readScalar(propsDict_.lookup("radius"))),
    tau_(readScalar(propsDict_.lookup("tau")))    
{

//     agentIds_.clear();
// 
//     selectAgentIds ids
//     (
//         cloud_.cP(),
//         propsDict_
//     );
// 
//     agentIds_ = ids.agentIds();
    
    initialTimeDelay_ = 0.0;
    
    if (propsDict_.found("initialTimeDelay"))
    {
        initialTimeDelay_ = readScalar(propsDict_.lookup("initialTimeDelay"));
    }
    
    initialTime_ = time_.timeOutputValue();
   
    
    // load tracking numbers
    
    trackingNumbers_ = List<List<label> >(propsDict_.lookup("trackingNumbers"));
    
    Info << "trackingNumbers = " << trackingNumbers_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentGroupForce::~agentGroupForce()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void agentGroupForce::initialConfiguration()
{
   
}

void agentGroupForce::force(agent* p2)
{
    if((time_.timeOutputValue() - initialTime_) > initialTimeDelay_)
    {
        // find group members 
        label tN2 = p2->trackingNumber();
        label groupId = -1;
        
        // being very inefficent here - any optimisation is welcome ( won't scale up well )
        forAll(trackingNumbers_, i)
        {
            label id = findIndex(trackingNumbers_[i], tN2);
            
            if(id != -1)
            {
               groupId = i;
            }
        }
            
        if(groupId != -1)
        {
            // find leader of group
            label tN1 = trackingNumbers_[groupId][0];
            
            if(tN1 != tN2)
            {
                if(cloud_.tracker().isTrackingNumberAvailable(tN1))
                {
                    agent* p1 = cloud_.tracker().getAgent();
                    vector r1 = p1->position();
                    
                    vector r2 = p2->position();
                    vector n = r1 - r2;
                    scalar magR12 = mag(n);
                    n /= magR12;
                    
                    if(magR12 > radius_)
                    {
                        p2->f() += (convergentSpeed_*n - p2->v())*p2->mass() / tau_; // MAKE SURE n is UNIT vector
                    }
                }
            }
        }
    }
}

void agentGroupForce::newForce()
{
    
}




} // End namespace Foam

// ************************************************************************* //
