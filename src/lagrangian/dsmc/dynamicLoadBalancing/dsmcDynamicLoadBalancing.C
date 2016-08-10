/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class
    dsmcDynamicLoadBalancing

Description

\*----------------------------------------------------------------------------*/

#include "dsmcDynamicLoadBalancing.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "wallPolyPatch.H"
#include "dsmcCloud.H"
#include <unistd.h>

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct for dsmcInitialise
// dsmcDynamicLoadBalancing::dsmcDynamicLoadBalancing
// (
// )
// :
// {}


//- Construct from mesh, cloud and boolean (dsmcFoam)
dsmcDynamicLoadBalancing::dsmcDynamicLoadBalancing
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    time_(t),
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    performBalance_(false),
    originalEndTime_(time_.time().endTime().value()),
    dsmcLoadBalanceDict_
    (
        IOobject
        (
            "loadBalanceDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    controlDict_
    (
        IOobject
        (
            "controlDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        )
    ),
    maxImbalance_(readScalar(dsmcLoadBalanceDict_.lookup("maximumAllowableImbalance")))
{}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDynamicLoadBalancing::~dsmcDynamicLoadBalancing()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDynamicLoadBalancing::clean()
{

}


void dsmcDynamicLoadBalancing::update
(
)
{
    if(time_.time().outputTime())
    {
        // Load Balancing
        
        Switch enableBalancing = dsmcLoadBalanceDict_.lookup("enableBalancing");
        
        if ( Pstream::parRun() )
        {
            const scalar& allowableImbalance = maxImbalance_;
                
            //First determine current level of imbalance - do this for all
            // parallel runs, even if balancing is disabled
            scalar nGlobalParticles = cloud_.size();
            Foam::reduce(nGlobalParticles, sumOp<scalar>());
            
            scalar idealNParticles = scalar(nGlobalParticles)/scalar(Pstream::nProcs());
            
            scalar nParticles = cloud_.size();
            scalar localImbalance = mag(nParticles - idealNParticles);
            Foam::reduce(localImbalance, maxOp<scalar>());
            scalar maxImbalance = localImbalance/idealNParticles;
            
            Info << "    Maximum imbalance = " << 100*maxImbalance << "%" << endl;
            
            if( maxImbalance > allowableImbalance && enableBalancing)
            {   
                performBalance_ = true;
                
                scalar currentTime = time_.time().value();

                if (Pstream::master())
                {                                       
                    controlDict_.set("endTime",currentTime);
                    controlDict_.Foam::regIOobject::write();
                                        
                    system("cp processor0/system/controlDict system/controlDict");
                }
            }
        }
    }
}
    
void dsmcDynamicLoadBalancing::perform
(
)
{
    if(performBalance_)
    {
        if (Pstream::master())
        {   
//             string userDir = system("echo $WM_PROJECT_USER_DIR");
//             string compileOptions1 = system("echo $WM_ARCH");
//             string compileOptions2 = system("echo $WM_COMPILER");
//             string compileOptions3 = system("echo $WM_COMPILE_OPTION");
//             
//             string reconstructPar = userDir+"/platforms/"+compileOptions1+compileOptions2+compileOptions3+"/bin/reconstructPar -latestTime";
            
//             Info << "reconstructPar = " << reconstructPar << endl;
            
            system("reconstructPar -latestTime");
                     
            system("decomposeDSMCLoadBalancePar -force");
//             system("/home/cwhite/OpenFOAM/OpenFOAM-2.4.0-MNF/platforms/linux64GccMNFDPOpt/bin/decomposeDSMCLoadBalancePar -force");
            performBalance_ = false;
            
            controlDict_.set("endTime",originalEndTime_);
            controlDict_.Foam::regIOobject::write();
                    
            system("cp processor0/system/controlDict system/controlDict");
        }
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

}  // End namespace Foam

// ************************************************************************* //
