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
    
    clearLagrangianFields();
    calculateForce();
    updateAcceleration();
    
    cloud_.fields().createFields();
    cloud_.controllers().initialConfig();
}

void velocityVerlet::evolve()
{
    cloud_.controllers().controlVelocitiesI();
    
    updateHalfVelocity();

    cloud_.controllers().controlBeforeMove();
    
    cloud_.move();
    
    cloud_.buildCellOccupancy();

    cloud_.b().afterMove();    

    cloud_.controllers().controlBeforeForces();
    
    clearLagrangianFields();
    
    calculateForce();
    
    updateAcceleration();
    
    cloud_.b().afterForce();  
    
    cloud_.controllers().controlAfterForces();
        
//     cloud_.controlAfterForces();
    updateHalfVelocity();
    
    cloud_.controllers().controlVelocitiesII(); 
    
//     test();
//     cloud_.controlAfterVelocity();
    
    postTimeStep();
}

void velocityVerlet::postTimeStep()
{
    cloud_.fields().calculateFields();
    cloud_.fields().writeFields();
    
    cloud_.controllers().calculateStateProps();
    cloud_.controllers().outputStateResults();    
}

// void velocityVerlet::test()
// {
//     IDLList<agent>::iterator mol(cloud_.begin());
//     
//     label i = 0;
//     
//     for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
//     {
//         if(i==10)
//         {
//             Info << "position = " << mol().position() 
//                  << ", velocity = " << mol().v()  
//                  << endl;
//         }
//         
//         i++;
//     }
// }

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

void velocityVerlet::calculateForce()
{
    cloud_.f().calculatePairForces();
}

void velocityVerlet::clearLagrangianFields()
{
    IDLList<agent>::iterator mol(cloud_.begin());

    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {
        mol().a() = vector::zero;

        mol().f() = vector::zero;

        mol().R() = GREAT;
    }
}

void velocityVerlet::updateAcceleration()
{
    IDLList<agent>::iterator mol(cloud_.begin());

    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {
        if(!mol().frozen())
        {
//             scalar mass = cloud_.cP().mass(mol().id());
                
            mol().a() = mol().f()/mol().mass();
        }
    }
}


} // End namespace Foam

// ************************************************************************* //
