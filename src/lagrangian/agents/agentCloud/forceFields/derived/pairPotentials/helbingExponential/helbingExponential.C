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
            
            thetaRef1_ = 90;
            thetaRef1_ *= constant::mathematical::pi/180;
            
            if (propsDict_.found("theta1"))
            {
                thetaRef1_ = readScalar(propsDict_.lookup("theta1"));
                thetaRef1_ *= constant::mathematical::pi/180;                
            }             
            
        }
//         else if(option == "opposingFlowAndAnisotropy")
//         {
//             option_ = 2;
//             
//             thetaRef2_=10; //deg
//             
//         
//             if (propsDict_.found("theta2"))
//             {
//                 thetaRef2_ = readScalar(propsDict_.lookup("theta2"));
//             }            
//             
//             C_ = readScalar(propsDict_.lookup("C"));
//             
//         }
    }
    
    injury_ = false;
    
    if (propsDict_.found("injury"))
    {
        injury_ = Switch(propsDict_.lookup("injury"));
        
        timeDelay_ = 0.0;
        
        if (propsDict_.found("injuryTimeDelay"))
        {
            timeDelay_ = readScalar(propsDict_.lookup("injuryTimeDelay"));
        }
        
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
    if(r > 0.0)
    {
        vector pairForce = vector::zero;
        vector socialForce = vector::zero;
        vector normalForce = vector::zero;
        vector frictionForce = vector::zero;
        
        vector rij = molI->position()-molJ->position();
        vector nij = rij/r;
        vector tij = vector(-nij.y(), nij.x(), 0);
        scalar dIJ = r;
        
        // the sum of radii
        scalar rIJ = molI->radius() + molJ->radius();
        
        if(dIJ > rIJ)
        {
            socialForce = A_*exp((rIJ-dIJ)/B_)*nij;
        } 
        else
        {
            socialForce = (A_*exp((rIJ-dIJ)/B_))*nij;
            normalForce = k_*(rIJ-dIJ)*nij;
            frictionForce = kappa_*(rIJ-dIJ)*((molJ->v() - molI->v()) & tij)*tij;
            
            if(injury_)
            {
                if(cloud_.mesh().time().timeOutputValue() > timeDelay_)
                {
                    if( (mag(normalForce)/(2*constant::mathematical::pi*molI->radius())) >= maxForce_)
                    {
                        molI->v() = vector::zero;
                        molI->special()= -2;
                        molI->id() = agentId_;
                    }
                    
                    if( (mag(normalForce)/(2*constant::mathematical::pi*molJ->radius())) >= maxForce_)
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
        }
        
        // apply force
        
        if(option_ == 0) // standard
        {
            pairForce = socialForce + normalForce  + frictionForce;            
            molI->f() += pairForce;
            molJ->f() += -pairForce;
        }
        else if(option_ == 1) //|| (option_ == 2)) // anisotropy
        {
            scalar wI = 0.0;

            scalar magVI = mag(molI->v());
            
            if(magVI > 0.0)
            {    
                vector vI = molI->v()/magVI;
                
                scalar dotI = vI & -nij;

                if(dotI > 0)
                {
                    if(dotI < 1)
                    {
                        wI = 1.0 - (acos(dotI)/thetaRef1_);
                    }
                    else
                    {
                        wI = 1.0;
                    }
                }
            }
            
            scalar wJ = 0.0;
            
            scalar magVJ = mag(molJ->v());
            
            if(magVJ)
            {
                vector vJ = molJ->v()/magVJ;
                
                scalar dotJ = vJ & nij;
                
                if(dotJ > 0)
                {
                    if(dotJ < 1)
                    {
                        wJ = 1.0 - (acos(dotJ)/thetaRef1_);
                    }
                    else
                    {
                        wJ = 1.0;
                    }
                }
            }
            
            molI->f() += wI*socialForce;
            molJ->f() += -wJ*socialForce;        
            
            molI->f() += normalForce  + frictionForce;
            molJ->f() += -(normalForce  + frictionForce);            
            /*
            if(option_ == 2) // opposing flow
            {
                scalar vIJ = molI->v() & molJ->v();
                wI =0.0;
                wJ =0.0;
                
                if(vIJ < 0)
                {
                    if(magVI > 0.0)
                    {
                        vector vI = molI->v()/magVI;
                        scalar dotI = vI & -nij;
                        
                        if( (dotI > 0) && (dotI < 1))
                        {
                            scalar theta = acos(dotI);
                            
                            wI = 1.0 - (theta*thetaRef1_*constant::mathematical::pi/180.0);
                        }
                    }
                    
                    if(magVJ > 0.0)
                    {
                        vector vJ = molJ->v()/magVJ;
                        scalar dotJ = vJ & nij;
                        
                        if( (dotJ > 0) && (dotJ < 1))
                        {
                            scalar theta = acos(dotJ);
                            
                            wJ = 1.0 - (theta*thetaRef1_*constant::mathematical::pi/180.0);
                        }
                    }                    
                }
                
                vector tI = vector(-molI->v().y(), molI->v().x(), 0);
                scalar magtI = mag(tI);
                
                if(magtI > 0.0)
                {
                    tI /= magtI;
                    
                    molI->f() += C_*wI*tI*mag(socialForce);

//                     Info << "pos I = " << molI->position()
//                         << ", vI = " << molI->v()
//                         << ", wI = " << wI 
//                         << ", socialForce = " << socialForce
//                         << ", tI = " << tI
//                         << ", rIJ = " << r
//                         << ", force = " << C_*wI*tI*mag(socialForce)
//                         << endl;
                    
                    
                    molJ->f() += -C_*wJ*tI*mag(socialForce);
                    
//                     Info << "pos J = " << molJ->position()
//                         << ", vJ = " << molJ->v()
//                         << ", wJ = " << wJ 
//                         << ", pairForce = " << socialForce
//                         << ", tJ = " << -tI
//                         << ", rIJ = " << r
//                         << ", force = " << -C_*wJ*tI*mag(socialForce)
//                         << endl;                     
                }
            }
            */
        }

    }
    else
    {
        Info << "WARNING: two agents are overlapping, so no force will be applied" << endl;
    }
    
}




} // End namespace Foam

// ************************************************************************* //
