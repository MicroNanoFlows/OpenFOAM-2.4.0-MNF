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

#include "scheduledForce.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(scheduledForce, 0);

addToRunTimeSelectionTable(bodyForce, scheduledForce, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
scheduledForce::scheduledForce
(
    agentCloud& cloud,
    Time& t,
    const dictionary& dict
)
:
    bodyForce(cloud, t, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    agentIds_()
{

    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();
    
    threshold_ = 0.1;
    
    // to be modified more elegantly in the future
    factor_ = readScalar(propsDict_.lookup("factor"));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

scheduledForce::~scheduledForce()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void scheduledForce::initialConfiguration()
{}

void scheduledForce::force(agent* p)
{
    if(findIndex(agentIds_, p->id()) != -1)
    {
        // if the destination is outside the bounds, there is no destination to go to
        if(p->d().x() >= 0.0)
        {
            vector rIJ = p->d() - p->position();
            
            scalar rIJMag = mag(rIJ);
            
            if(rIJMag < threshold_)
            {
                //arrived at its destination - reset it to outside domain
                p->d() = vector(-1,-1,-1);
            }
            else
            {
                // ideally we need to insert route based selection here
                vector n = rIJ/rIJMag; // unit vector
                
                p->f() += n*factor_;
            }
        }
    }
}

void scheduledForce::newForce()
{
    
}





} // End namespace Foam

// ************************************************************************* //
