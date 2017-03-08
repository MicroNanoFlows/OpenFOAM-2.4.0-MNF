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

#include "debuggerAgent.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(debuggerAgent, 0);

addToRunTimeSelectionTable(agentController, debuggerAgent, dictionary);





// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
debuggerAgent::debuggerAgent
(
    Time& t,
    agentCloud& cloud,
    const dictionary& dict
)
:
    agentController(t,  cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
//     model_(),
    agentIds_(),
    bb_(mesh_.bounds().min(), mesh_.bounds().max()),
    nTimeSteps_(0.0),
    force_(vector::zero)
   
{
//     writeInTimeDir_ = true;
//     writeInCase_ = true;


    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();

    trackingNumbers_ = List<label>(propsDict_.lookup("trackingNumbers"));
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

debuggerAgent::~debuggerAgent()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void debuggerAgent::initialConfiguration()
{}

void debuggerAgent::controlBeforeVelocityI()
{}

void debuggerAgent::controlBeforeMove()
{}

void debuggerAgent::controlBeforeForces()
{}

void debuggerAgent::controlDuringForces
(
    agent* molI,
    agent* molJ
)
{}

void debuggerAgent::controlAfterForces()
{
//     Info << "debuggerAgent"  << endl;
// 
//     IDLList<agent>::iterator mol(cloud_.begin());
// 
//     for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
//     {
//         if(findIndex(agentIds_, mol().id()) != -1)
//         {
//             if(bb_.contains(mol().position()))
//             {
//                 
//             }
//             else
//             {
//                 FatalErrorIn("debuggerAgent") << nl
//                     << " Agent position just left the domain at : "
//                     << mol().position() 
//                     << ", tracking Number = " << mol().trackingNumber()
//                     << nl << abort(FatalError);                         
//             }
//         }
//     }
}


void debuggerAgent::controlAfterVelocityII()
{}

void debuggerAgent::calculateProperties()
{
    Info << "debuggerAgent"  << endl;

    IDLList<agent>::iterator mol(cloud_.begin());

    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {
        if(findIndex(trackingNumbers_, mol().trackingNumber()) != -1)
        {
//             if(bb_.contains(mol().position()))
//             {
//                 
//             }
//             else
//             {
//                 FatalErrorIn("debuggerAgent") << nl
//                     << " Agent position just left the domain at : "
//                     << mol().position() 
//                     << ", tracking Number = " << mol().trackingNumber()
//                     << nl << abort(FatalError);                         
//             }
            scalar Vmag = mag(mol().v());
            
            Info << "tN = " << mol().trackingNumber() 
                << ", pos = " << mol().position()
                << ", F = " << mol().f()
                << ", a = " << mol().a()
                << ", v = " << mol().v()
                << ", vMag = " << Vmag
                << ", desired speed = " << mol().desiredSpeed()
                << ", destination = " << mol().d()
                << ", time = " << mol().t()
                << endl;
            
        }
    }   
    
    
}

void debuggerAgent::output
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
