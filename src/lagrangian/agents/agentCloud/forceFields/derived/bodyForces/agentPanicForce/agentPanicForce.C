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

#include "agentPanicForce.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agentPanicForce, 0);

addToRunTimeSelectionTable(bodyForce, agentPanicForce, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
agentPanicForce::agentPanicForce
(
    agentCloud& cloud,
    Time& t,
    const dictionary& dict
)
:
    bodyForce(cloud, t, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    agentIds_(),
//     desiredSpeed_(readScalar(propsDict_.lookup("desiredSpeed"))),
//     desiredDirection_(propsDict_.lookup("desiredDirection")),
    tau_(readScalar(propsDict_.lookup("tau")))
    
{

    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();
    
    initialTimeDelay_ = 0.0;
    
    if (propsDict_.found("initialTimeDelay"))
    {
        initialTimeDelay_ = readScalar(propsDict_.lookup("initialTimeDelay"));
    }
    
    initialTime_ = time_.timeOutputValue();
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentPanicForce::~agentPanicForce()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void agentPanicForce::initialConfiguration()
{
   
}

void agentPanicForce::force(agent* p)
{
    if((time_.timeOutputValue() - initialTime_) > initialTimeDelay_)
    {
        if(findIndex(agentIds_, p->id()) != -1)
        {
            if(mag(p->f()) < 0.01)
            {  
                p->f() += (p->desiredSpeed()*n - p->v())*p->mass() / tau_; // MAKE SURE n is UNIT vector
            }
        }
    }
}

void agentPanicForce::newForce()
{
    
}






} // End namespace Foam

// ************************************************************************* //
