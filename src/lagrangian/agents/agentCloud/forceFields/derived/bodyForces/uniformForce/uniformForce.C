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

#include "uniformForce.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(uniformForce, 0);

addToRunTimeSelectionTable(bodyForce, uniformForce, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
uniformForce::uniformForce
(
    agentCloud& cloud,
    Time& t,
    const dictionary& dict
)
:
    bodyForce(cloud, t, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    model_(),
    agentIds_()
{
    model_ = autoPtr<gravityForce>
    (
        gravityForce::New(t, propsDict_)
    );

    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();

    setBoundBoxes();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

uniformForce::~uniformForce()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void uniformForce::initialConfiguration()
{}

void uniformForce::force(agent* p)
{
    if(findIndex(agentIds_, p->id()) != -1)
    {
        forAll(boxes_, b)
        {
            if(boxes_[b].contains(p->position()))
            {
                vector force = vector::zero;
                
                if(model_->timeVarying())   
                {
                    const scalar t = time_.timeOutputValue();
                    
                    force = model_->force(t);
                }
                else if(model_->spaceVarying())
                {
                    force = model_->force(p->position());
                }
                
                p->f() += force;
                
                p->a() += force/p->mass();
            }
        }
    }
}

void uniformForce::newForce()
{
    
}


void uniformForce::setBoundBoxes()
{
 
    PtrList<entry> boxList(propsDict_.lookup("boxes"));

    boxes_.setSize(boxList.size());

    forAll(boxList, b)
    {
        const entry& boxI = boxList[b];
        const dictionary& dict = boxI.dict();

        vector startPoint = dict.lookup("startPoint");
        vector endPoint = dict.lookup("endPoint");
        boxes_[b].resetBoundedBox(startPoint, endPoint);
    }
}



} // End namespace Foam

// ************************************************************************* //
