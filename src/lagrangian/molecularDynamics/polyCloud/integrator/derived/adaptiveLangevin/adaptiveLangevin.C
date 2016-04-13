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

#include "adaptiveLangevin.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(adaptiveLangevin, 0);

addToRunTimeSelectionTable(polyIntegrator, adaptiveLangevin, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
adaptiveLangevin::adaptiveLangevin
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyIntegrator(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    gamma_(readScalar(propsDict_.lookup("gamma"))),
    T_(readScalar(propsDict_.lookup("temperature")))    
{
    deltaT_ = mesh_.time().deltaT().value();
    xi_ = gamma_;
    
    const scalar& kB = molCloud_.redUnits().kB();
    kB_ = kB;
    
    sigma_ = sqrt(2.0*gamma_*kB_*T_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

adaptiveLangevin::~adaptiveLangevin()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void adaptiveLangevin::init()
{
    
    scalar dof = 0.0;
    scalar nMols = 0.0;
    
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(!mol().frozen())
            {
                dof += molCloud_.cP().degreesOfFreedom(mol().id());
                nMols += 1.0;
            }
        }    
    }   
    
    if(Pstream::parRun())
    {
        reduce(nMols, sumOp<scalar>());        
        reduce(dof, sumOp<scalar>());
    }
        
    Nd_ = dof;
    N_ = nMols;
    
    scalar mu = N_; 
    invMu_= 1/mu;    
}


void adaptiveLangevin::evolve()
{
    updateHalfVelocity();
    
    {
        scalar kE = 0.0;
        scalar angularKeSum = 0.0;
        calculateKE(kE, angularKeSum);
        kE_.append(kE);
        angularKE_.append(angularKeSum);        
    }    
    
    // move 1
    molCloud_.move(0.5*deltaT_);
    molCloud_.updateAfterMove(0.5*deltaT_);
    
    molCloud_.buildCellOccupancy();
    
    calculateXi();
    
    // update velocity 
    if(xi_ != 0)
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(!mol().frozen())
            {
                const scalar& mass = molCloud_.cP().mass(mol().id());
                
                vector rand = vector::zero;
                
                rand.x() = molCloud_.rndGen().sample01<scalar>();
                rand.y() = molCloud_.rndGen().sample01<scalar>();
                rand.z() = molCloud_.rndGen().sample01<scalar>();
                
                // double check
                scalar e=exp(-xi_*deltaT_);
                mol().v() = (e*mol().v()) + sigma_*sqrt((1-e*e)/(2*xi_*mass))*rand;
            }
        }    
    }
    else
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(!mol().frozen())
            {
                const scalar& mass = molCloud_.cP().mass(mol().id());

                vector rand = vector::zero;
                
                rand.x() = molCloud_.rndGen().sample01<scalar>();
                rand.y() = molCloud_.rndGen().sample01<scalar>();
                rand.z() = molCloud_.rndGen().sample01<scalar>();

                mol().v() += sqrt(deltaT_/mass)*sigma_*rand;
            }
        }        
    }
    
    // calculate  KE
    {
        scalar kE = 0.0;
        scalar angularKeSum = 0.0;
        calculateKE(kE, angularKeSum);
        kE_.append(kE);
        angularKE_.append(angularKeSum);        
    }
        
    calculateXi();
    
    // move 2
    molCloud_.move(0.5*deltaT_);
    molCloud_.updateAfterMove(0.5*deltaT_);
    
    
    molCloud_.buildCellOccupancy();
    molCloud_.clearLagrangianFields();
    molCloud_.calculateForce();
    molCloud_.updateAcceleration();
        
    updateHalfVelocity();
}

void adaptiveLangevin::updateHalfVelocity()
{
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(!mol().frozen())
        {
            mol().v() += 0.5*mol().a()*deltaT_;
        }
    }
}

void adaptiveLangevin::calculateXi()
{
    scalar Vtot = 0.0;
    
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(!mol().frozen())
            {
                const scalar& mass = molCloud_.cP().mass(mol().id());
                
                Vtot += mass*(mol().v().x()*mol().v().x()) +
                        (mol().v().y()*mol().v().y()) +
                        (mol().v().z()*mol().v().z());
            }
        }
    }
    
    
    if(Pstream::parRun())
    {
        reduce(Vtot, sumOp<scalar>());        
    }        
    
    xi_ += 0.5*deltaT_*invMu_*(Vtot - Nd_*kB_*T_);    
}

void adaptiveLangevin::calculateKE(scalar& kE, scalar& angularKeSum)
{
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(!mol().frozen())
        {
            const scalar& mass = molCloud_.cP().mass(mol().id());
            
            kE += mass*magSqr(mol().v());

            const diagTensor& molMoI(molCloud_.cP().momentOfInertia(mol().id()));

            // angular speed 
            const vector& molOmega(inv(molMoI) & mol().pi());                
            
            angularKeSum += 0.5*(molOmega & molMoI & molOmega);
        }
    }    
    
    if(Pstream::parRun())
    {
        reduce(kE, sumOp<scalar>());        
        reduce(angularKeSum, sumOp<scalar>());
    }
}


} // End namespace Foam

// ************************************************************************* //
