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

#include "helbingVelocityDependentBAnisotropy.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(helbingVelocityDependentBAnisotropy, 0);

addToRunTimeSelectionTable(pairPotential, helbingVelocityDependentBAnisotropy, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
helbingVelocityDependentBAnisotropy::helbingVelocityDependentBAnisotropy
(
    agentCloud& cloud,
    const word& name,
    const dictionary& dict
)
:
    pairPotential(cloud, name, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    A_(readScalar(propsDict_.lookup("A"))),
    Bmax_(readScalar(propsDict_.lookup("Bmax"))),
    Bmin_(readScalar(propsDict_.lookup("Bmin"))),
    k_(readScalar(propsDict_.lookup("k"))),
    kappa_(readScalar(propsDict_.lookup("kappa"))),
    lambda_(readScalar(propsDict_.lookup("lambda")))
{

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

helbingVelocityDependentBAnisotropy::~helbingVelocityDependentBAnisotropy()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void helbingVelocityDependentBAnisotropy::initialConfiguration()
{}

scalar helbingVelocityDependentBAnisotropy::energy(const scalar& r)
{
    scalar energy = 0.0;

    return energy;
}
        
scalar helbingVelocityDependentBAnisotropy::force(const scalar& r)
{
    scalar force = 0.0;
    
    return force;
}

void helbingVelocityDependentBAnisotropy::pairPotentialFunction
(
    agent* molI,
    agent* molJ,
    const scalar& r,
    scalar& energy,
    vector& force
)
{
    // r is the distance between two agents 
    
    vector rij = molI->position()-molJ->position();
//     scalar rijMag=mag(rij);
    vector nij = rij/r;
    vector tij = vector (-nij.y(), nij.x(), 0);
    
    scalar dIJ = r;
    
    // the sum of radii
    scalar rIJ = molI->radius() + molJ->radius();
    
    // Read Bmax and Bmin and define linear dependence w.r.t. velocity
    scalar Bvel = Bmin_ + (Bmax_ - Bmin_)/1.5 * mag(molI->v());
    
    
    //Add anisotropy parameter!
//     vector vAlpha = molI->v()/mag(molI->v());
//     scalar test = mag(vAlpha);
//     cout << ">>>>>> test is equal to: " << test; 
// 
//     vector dAlphaBeta = (-1)*rij/mag(rij);
//     scalar cosAngle = vAlpha & dAlphaBeta;
//     cout << ">>>>>> cosAngle is equal to: " << cosAngle; 
//     scalar w = (lambda_ + (1 - lambda_) * (1+cosAngle)/2);
//     cout << ">>>>>> w is equal to: " << w; 
    
    scalar w;
    if(mag(molI->v()) == 0)
    {
	w = 0;
    }
    else
    {
	vector v1 = molI->v()/mag(molI->v());
	vector d12 = (-1)*rij/mag(rij);
	scalar cosAngle = v1 & d12;
	w = (lambda_ + (1 - lambda_) * (1+cosAngle)/2);
    }
    
    
    if(dIJ >= rIJ)
    {
        force = w*A_*exp((rIJ-dIJ)/Bvel)*nij;
    } 
    else
    {
        force = (w*A_*exp((rIJ-dIJ)/Bvel) + k_*(rIJ-dIJ) )*nij  + 
                kappa_*(rIJ-dIJ)*((molJ->v() - molI->v()) & tij)*tij;
    }
}




} // End namespace Foam

// ************************************************************************* //
