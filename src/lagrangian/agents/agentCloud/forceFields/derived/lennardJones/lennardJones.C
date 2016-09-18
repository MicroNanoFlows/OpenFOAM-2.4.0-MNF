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

#include "lennardJones.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(lennardJones, 0);

addToRunTimeSelectionTable(pairPotential, lennardJones, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
lennardJones::lennardJones
(
    agentCloud& cloud,
    const word& name,
    const dictionary& dict
)
:
    pairPotential(cloud, name, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    sigma_(readScalar(propsDict_.lookup("sigma"))),
    epsilon_(readScalar(propsDict_.lookup("epsilon")))      
{
    if(rU_.runReducedUnits())
    {
        sigma_ /= rU_.refLength();
        epsilon_ /= rU_.refEnergy();
    }
    
    E_at_Rmin_ = energy(rMin_);        
    F_at_Rmin_ = force(rMin_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

lennardJones::~lennardJones()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void lennardJones::initialConfiguration()
{}

scalar lennardJones::energy(const scalar& r)
{
    // (rIJ/sigma)^-2
    scalar ir2 = (sigma_/r)*(sigma_/r);

    // (rIJ/sigma)^-6
    scalar ir6 = ir2*ir2*ir2;

    return 4.0*epsilon_*(ir6*(ir6 - 1.0));
}
        
scalar lennardJones::force(const scalar& r)
{
    // (rIJ/sigma)^-2
    scalar ir2 = (sigma_/r)*(sigma_/r);

    // (rIJ/sigma)^-6
    scalar ir6 = ir2*ir2*ir2;    
        
    return 24.0*epsilon_*ir6*(2.0*ir6 - 1.0)*sqrt(ir2);
}

void lennardJones::pairPotentialFunction
(
    agent* molI,
    agent* molJ,
    const scalar& r,
    scalar& energy,
    scalar& force
)
{
    // (rIJ/sigma)^-2
    scalar ir2 = (sigma_/r)*(sigma_/r);

    // (rIJ/sigma)^-6
    scalar ir6 = ir2*ir2*ir2;  
    
    energy = E_at_Rmin_;    
    force = F_at_Rmin_;

    
    if(r > rMin_)
    {
        energy = 4.0*epsilon_*(ir6*(ir6 - 1.0));    
        force = 24.0*epsilon_*ir6*(2.0*ir6 - 1.0)*sqrt(ir2);            
    }
}




} // End namespace Foam

// ************************************************************************* //
