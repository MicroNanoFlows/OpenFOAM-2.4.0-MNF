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

#include "agentWillForceScheduled.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agentWillForceScheduled, 0);

addToRunTimeSelectionTable(bodyForce, agentWillForceScheduled, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
agentWillForceScheduled::agentWillForceScheduled
(
    agentCloud& cloud,
    Time& t,
    const dictionary& dict
)
:
    bodyForce(cloud, t, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    agentIds_(),
//     desiredSpeed_(readScalar(propsDict_.lookup("desiredSpeed"))),
//     desiredDirection_(propsDict_.lookup("desiredDirection")),
    tau_(readScalar(propsDict_.lookup("tau"))),
    deltaT_(time_.deltaT().value())
    
{

    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();
    
    initialTimeDelay_ = 0.0;
    
    if (propsDict_.found("initialTimeDelay"))
    {
        initialTimeDelay_ = readScalar(propsDict_.lookup("initialTimeDelay"));
    }
    else
    {
	initialTimeDelay_ = 0;
    }
    
/*    panic_ = false;
    
    if (propsDict_.found("panic"))
    {
        panic_ = Switch(propsDict_.lookup("panic"));
        smallForce_= readScalar(propsDict_.lookup("smallForce"));
        forceMag_= readScalar(propsDict_.lookup("forceMag"));
        
    } */   
    
    initialTime_ = time_.timeOutputValue();
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentWillForceScheduled::~agentWillForceScheduled()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void agentWillForceScheduled::initialConfiguration()
{
   
}

void agentWillForceScheduled::force(agent* p)
{
        
//     if((time_.timeOutputValue() - initialTime_) > initialTimeDelay_)
//     {
//         if(findIndex(agentIds_, p->id()) != -1)
//         {
//             
//             if(p->t() == 0.0)
//             {  
//                 vector n = p->d() - p->position();
//                 n /= mag(n);
//                 
//                 p->f() += (p->desiredSpeed()*n - p->v())*p->mass() / tau_; // MAKE SURE n is UNIT vector
//                 
// //                 if(panic_)
// //                 {
// //                     if(mag(p->f()) < smallForce_)
// //                     {
// //                         Info << "force = " << mag(p->f()) << endl;
// //                         
// //                         p->f() = n*forceMag_;
// //                     }
// //                 }
// //                 if(p->trackingNumber() == 75 )
// //                 {
// //                     Info << "75, force = " <<  p->f() << endl;
// //                 }
//             }
//         }
//     }
//     else
//     {
//         p->v() = vector::zero;
//     }
    
    
//     Info << "deltaT_= " << deltaT_ << endl;
    
    if((time_.timeOutputValue() - initialTime_) > initialTimeDelay_)
    {
        if(findIndex(agentIds_, p->id()) != -1)
        {
            
            if(p->t() <= 0.0)
            {  
                vector n = p->d() - p->position();
                n /= mag(n);
                
                p->f() += (p->desiredSpeed()*n - p->v())*p->mass() / tau_;
            }
            else
	    {
// 		p->v() = vector::zero;
// 		p->t() -=deltaT_;
// 		Info << "p->(t) = " << p->t() << endl;
	    }
        }
    }
    else
    {
//         p->v() = vector::zero;
	
	if ((time_.timeOutputValue() - initialTime_) <= 0 && initialTimeDelay_ > 0)
	{
	    scalar uniformTimeDelay = initialTimeDelay_ + (30-initialTimeDelay_)*cloud_.rndGen().scalar01();
	    p->t() = uniformTimeDelay;
	    Info << "p->(t) = " << p->t() << endl;
	}
	
// 	p->t() -= 0.5*deltaT_;
// 	Info << "p->(t) = " << p->t() << endl;
	
    }
    
    
    
    
    
}

void agentWillForceScheduled::newForce()
{
    
}






} // End namespace Foam

// ************************************************************************* //
