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

#include "velocityVerlet.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(velocityVerlet, 0);

addToRunTimeSelectionTable(agentIntegrator, velocityVerlet, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
velocityVerlet::velocityVerlet
(
    Time& t,
    agentCloud& cloud,
    const dictionary& dict
)
:
    agentIntegrator(t, cloud, dict)
/*    propsDict_(dict.subDict(typeName + "Properties")),*/    
  
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

velocityVerlet::~velocityVerlet()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void velocityVerlet::init()
{
    cloud_.f().initialConfig();
    cloud_.b().initialConfig();
    cloud_.ob().initialConfig();        
    cloud_.s().initialConfig();
    
    cloud_.tracker().track();
    
    // initial force
    clearLagrangianFields();
    calculateForce();
    
    cloud_.fields().createFields();
    cloud_.controllers().initialConfig();
}

void velocityVerlet::evolve()
{
//     setDesiredDirection();
    
    cloud_.s().setSchedules();
    
    cloud_.controllers().controlVelocitiesI();
    
    updateHalfVelocity();
    
    cloud_.controllers().controlBeforeMove();

    checkMaxVelocity();
    
    cloud_.move();

    cloud_.tracker().track();
    
    cloud_.buildCellOccupancy();

    cloud_.b().afterMove();

    cloud_.controllers().controlBeforeForces();
    
    clearLagrangianFields();

    calculateForce();
    
    cloud_.controllers().controlAfterForces();
        
    updateHalfVelocity();
    
    cloud_.controllers().controlVelocitiesII(); 

    // after time step 
    
    cloud_.fields().calculateFields();
    cloud_.fields().writeFields();
    cloud_.controllers().calculateStateProps();
    cloud_.controllers().outputStateResults(); 
    cloud_.ob().calculateProperties();
    cloud_.ob().outputResults();
}






void velocityVerlet::updateHalfVelocity()
{
    scalar trackTime = mesh_.time().deltaT().value();
        
    IDLList<agent>::iterator mol(cloud_.begin());

    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {
        if(!mol().frozen())
        {
            mol().v() += 0.5*trackTime*mol().a();
        }
    }
}

void velocityVerlet::clearLagrangianFields()
{
    IDLList<agent>::iterator mol(cloud_.begin());

    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {
        mol().f() = vector::zero;

        mol().R() = GREAT;
    }
}

void velocityVerlet::calculateForce()
{
    cloud_.f().calculatePairForces(); // social forces
    
    cloud_.f().calculateBodyForces(); // body forces
    
    cloud_.b().afterForce(); // wall forces
    cloud_.f().calculateWallForces();
}

void velocityVerlet::checkMaxVelocity()
{
    IDLList<agent>::iterator mol(cloud_.begin());

    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {
        scalar vMax = cloud_.cP().vMax()[mol().id()];
        
        scalar magV = mag(mol().v());
        
        if(magV > vMax)
        {
            vector n = mol().v()/magV;
            mol().v() = vMax*n;
        }
    }
}




} // End namespace Foam

// ************************************************************************* //
