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

#include "helbingExponential.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(helbingExponential, 0);

addToRunTimeSelectionTable(pairPotential, helbingExponential, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
helbingExponential::helbingExponential
(
    agentCloud& cloud,
    const word& name,
    const dictionary& dict
)
:
    pairPotential(cloud, name, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    A_(readScalar(propsDict_.lookup("A"))),
    B_(readScalar(propsDict_.lookup("B"))),
    k_(readScalar(propsDict_.lookup("k"))),
    kappa_(readScalar(propsDict_.lookup("kappa")))
{

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

helbingExponential::~helbingExponential()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void helbingExponential::initialConfiguration()
{}

scalar helbingExponential::energy(const scalar& r)
{
    scalar energy = 0.0;

    return energy;
}
        
scalar helbingExponential::force(const scalar& r)
{
    scalar force = 0.0;
    
    return force;
}

void helbingExponential::pairPotentialFunction
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
    
    if(dIJ >= rIJ)
    {
        force = A_*exp((rIJ-dIJ)/B_)*nij;
    } 
    else
    {
        force = (A_*exp((rIJ-dIJ)/B_) + k_*(rIJ-dIJ) )*nij  + 
                kappa_*(rIJ-dIJ)*((molJ->v() - molI->v()) & tij)*tij;
    }
}




} // End namespace Foam

// ************************************************************************* //
