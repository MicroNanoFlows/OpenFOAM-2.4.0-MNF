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
    option_ = 0;
    
    if (propsDict_.found("option"))
    {
        const word option = propsDict_.lookup("option");
        
        if(option == "default")
        {
            option_ = 0;
        }
        else if(option == "anisotropic")
        {
            option_ = 1;
        }
    }
    
    injury_ = false;
    
    if (propsDict_.found("injury"))
    {
        injury_ = Switch(propsDict_.lookup("injury"));
        
        const word idName(propsDict_.lookup("agentIdInjured")); 
        const List<word>& idList(cloud_.cP().agentIds());

        label id = findIndex(idList, idName);

        if(id == -1)
        {
            FatalErrorIn("boxInitialise::setInitialConfiguration()")
                << "Cannot find molecule id: " << idName << nl << "in idList."
                << exit(FatalError);
        }        
        
        agentId_ = id;
        
        maxForce_ = readScalar(propsDict_.lookup("maxForce"));
    }

    
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
    
    if(dIJ > rIJ)
    {
        force = A_*exp((rIJ-dIJ)/B_)*nij;
    } 
    else
    {
        vector socialForce = (A_*exp((rIJ-dIJ)/B_))*nij;
        vector normalForce = k_*(rIJ-dIJ)*nij;
        vector frictionForce = kappa_*(rIJ-dIJ)*((molJ->v() - molI->v()) & tij)*tij;
        force = socialForce + normalForce  + frictionForce;
        
        if(injury_)
        {
            if( (mag(normalForce)*2*constant::mathematical::pi*molI->radius()) >= maxForce_)
            {
/*                Info << "position = " << molI->position()
                     <<  ", trackingNumber = " << molI->trackingNumber()
                    << " radius I = " << molI->radius()
                    <<  ", nij " << nij
                    <<  ", mag(nij) " << mag(nij)
                    <<  ", dIJ " << dIJ
                    <<  ", rIJ " << rIJ
                     << ", normalForce = " << mag(normalForce)
                     << ", force per unit width = " << (mag(normalForce)*2*constant::mathematical::pi*molI->radius())
                     << endl; */               

                     
                molI->v() = vector::zero;
                molI->special()= -2;
                molI->id() = agentId_;
            }
            
            if( (mag(normalForce)*2*constant::mathematical::pi*molJ->radius()) >= maxForce_)
            {
/*                Info << "position = " << molJ->position()
                    <<  ", trackingNumber = " << molJ->trackingNumber()
                    << " radius J = " << molJ->radius()
                    <<  ", nij " << nij
                    <<  ", mag(nij) " << mag(nij)
                    <<  ", dIJ " << dIJ
                    <<  ", rIJ " << rIJ                    
                     << ", normalForce = " << mag(normalForce)
                     << ", force per unit width = " << (mag(normalForce)*2*constant::mathematical::pi*molJ->radius())
                     << endl;*/                
                
                molJ->v() = vector::zero;
                molJ->special()= -2;
                molJ->id() = agentId_;
            }            
        }
                
    }
    
/*    scalar f = molI->fraction();

    if(molJ->fraction() < f)
    {
        f = molJ->fraction();
    } */   
    
    // apply force
    
    if(option_ == 0) // standard
    {
        molI->f() += force;
        molJ->f() += -force;
    }
    else if(option_ == 1) // ansitropy 
    {
        scalar wI = 0.0;
        scalar rDI = molI->dir() & -rij;
        
        if(rDI > 0)
        {
            scalar theta = acos(rDI/(mag(molI->dir())*mag(-rij)));
            
            wI = 1.0 - (2.0*theta/constant::mathematical::pi);
        }

        scalar wJ = 0.0;
        scalar rDJ = molJ->dir() & rij;
        
        if(rDJ > 0)
        {
            scalar theta = acos(rDJ/(mag(molJ->dir())*mag(rij)));
            
            wJ = 1.0 - (2.0*theta/constant::mathematical::pi);
        }
        
        molI->f() += wI*force;
        molJ->f() += -wJ*force;        
    }
}




} // End namespace Foam

// ************************************************************************* //
