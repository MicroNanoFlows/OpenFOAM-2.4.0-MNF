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

#include "agentWillForce.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agentWillForce, 0);

addToRunTimeSelectionTable(bodyForce, agentWillForce, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
agentWillForce::agentWillForce
(
    agentCloud& cloud,
    Time& t,
    const dictionary& dict
)
:
    bodyForce(cloud, t, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    agentIds_(),
    desiredSpeed_(readScalar(propsDict_.lookup("desiredSpeed"))),
    desiredDirection_(propsDict_.lookup("desiredDirection")),
    tau_(readScalar(propsDict_.lookup("tau")))
//     stdev_(readScalar(propsDict_.lookup("stdev")))
{

    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentWillForce::~agentWillForce()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void agentWillForce::initialConfiguration()
{
   
}

void agentWillForce::force(agent* p)
{
    if(findIndex(agentIds_, p->id()) != -1)
    {
        p->f() += (desiredSpeed_*desiredDirection_ - p->v())*p->mass() / tau_;
    }
}

void agentWillForce::newForce()
{
    
}




} // End namespace Foam

// ************************************************************************* //
