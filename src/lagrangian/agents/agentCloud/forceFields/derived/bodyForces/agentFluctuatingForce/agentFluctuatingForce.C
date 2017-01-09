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

#include "agentFluctuatingForce.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agentFluctuatingForce, 0);

addToRunTimeSelectionTable(bodyForce, agentFluctuatingForce, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
agentFluctuatingForce::agentFluctuatingForce
(
    agentCloud& cloud,
    Time& t,
    const dictionary& dict
)
:
    bodyForce(cloud, t, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    agentIds_(),
    gamma_(readScalar(propsDict_.lookup("gamma"))),
    T_(readScalar(propsDict_.lookup("temperature")))
{

    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();
    
    deltaT_ = time_.time().deltaT().value();
    
    kB_ = 1.38064852e-23;
    
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentFluctuatingForce::~agentFluctuatingForce()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void agentFluctuatingForce::initialConfiguration()
{}

void agentFluctuatingForce::force(agent* p)
{
    if(findIndex(agentIds_, p->id()) != -1)
    {
        scalar massI = p->mass();
        
        scalar factor = sqrt(2*gamma_*massI*kB_*T_/deltaT_);
        
        vector f = factor*vector
        (
            cloud_.rndGen().GaussNormal(),
            cloud_.rndGen().GaussNormal(),
            0.0
        ) - gamma_*massI*p->v();

//         Info << "factor = " << factor 
//             << ", force = " << f << endl;
        
        p->f() += f;
    }
}

void agentFluctuatingForce::newForce()
{
    
}




} // End namespace Foam

// ************************************************************************* //
