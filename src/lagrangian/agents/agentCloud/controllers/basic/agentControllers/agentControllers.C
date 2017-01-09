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

#include "agentControllers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//- Null Constructor 
agentControllers::agentControllers
(
    Time& t,
    const polyMesh& mesh
)
:
    time_(t),
    mesh_(mesh),
    agentControllersDict_
    (
        IOobject
        (
            "controllersDict",
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    nStateControllers_(0),
    stateControllersList_(),
    sCNames_(),
    sCIds_(),
    sCFixedPathNames_(),
    stateControllers_(),

    controllersDuringForceComp_()
{}


agentControllers::agentControllers
(
    Time& t,
    const polyMesh& mesh,
    agentCloud& cloud
)
:
    time_(t),
    mesh_(mesh),
    agentControllersDict_
    (
        IOobject
        (
            "controllersDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nStateControllers_(0),
    stateControllersList_(agentControllersDict_.lookup("controllers")),
    sCNames_(stateControllersList_.size()),
    sCIds_(stateControllersList_.size()),
    sCFixedPathNames_(stateControllersList_.size()),
    stateControllers_(stateControllersList_.size()),

    controllersDuringForceComp_()
{

    Info << nl << "Creating agentControllers" << nl << endl;

    //- state agentControllers

    DynamicList<label> controllersDuringForceComp(0);

    if(stateControllers_.size() > 0 )
    {
        forAll(stateControllers_, sC)
        {
            const entry& agentControllersI = stateControllersList_[sC];
            const dictionary& agentControllersIDict = agentControllersI.dict();
    
            stateControllers_[sC] = autoPtr<agentController>
            (
                agentController::New(time_, cloud, agentControllersIDict)
            );
    
            sCNames_[sC] = stateControllers_[sC]->type();
            sCIds_[sC] = sC;
    
            nStateControllers_++;

            if(stateControllers_[sC]->controlInterForces())
            {
                controllersDuringForceComp.append(sC);
            }
        }
    }

    //controllersDuringForceComp_.transfer(controllersDuringForceComp.shrink());
    controllersDuringForceComp_.transfer(controllersDuringForceComp);


    // creating directories for state controllers
    if(nStateControllers_ > 0)
    {
        // directory: case/controllers
        fileName controllersPath(time_.path()/"controllers");

        if( !isDir(controllersPath) )
        {
            mkDir(controllersPath);
        }
        else
        {
            rmDir(controllersPath);
            mkDir(controllersPath);             
        }

        forAll(stateControllers_, sC)
        {
            if(stateControllers_[sC]->writeInCase())
            {
                // directory: case/controllers/poly/stateControllers/<stateControllerModel>
                fileName stateControllerPath(controllersPath/sCNames_[sC]);

                if (!isDir(stateControllerPath))
                {
                    mkDir(stateControllerPath);    
                }
    
                const word& regionName = stateControllers_[sC]->regionName();

                // directory: case/controllers/poly/stateControllers/<stateControllerModel>/<cellZoneName>    
                fileName zonePath(stateControllerPath/regionName);
   
                if (!isDir(zonePath))
                {
                    mkDir(zonePath);    
                }
    
                sCFixedPathNames_[sC] = zonePath;
            }
        }
    }
}

agentControllers::~agentControllers()
{}

//- initial configuration
//- call this function after the agentCloud is completely initialised
void agentControllers::initialConfig()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->initialConfiguration();
    }
}

void agentControllers::controlVelocitiesI()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlBeforeVelocityI();
    }
}

void agentControllers::controlBeforeMove()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlBeforeMove();
    }    
}


void agentControllers::controlBeforeForces()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlBeforeForces();
    }
}


void agentControllers::controlDuringForceComputation
(
    agent* molI,
    agent* molJ
)
{
    forAll(controllersDuringForceComp_, n)
    {
        const label& sC = controllersDuringForceComp_[n];
        stateControllers_[sC]->controlDuringForces(molI, molJ);
    }
}

//- control molecular state -- call this after the intermolecular force calulation
void agentControllers::controlAfterForces()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlAfterForces();
    }
}

void agentControllers::controlVelocitiesII()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlAfterVelocityII();
    }
}

//- calculate properties -- call this at the end of the MD time-step.
void agentControllers::calculateStateProps()
{
    forAll(stateControllers_, sC)
    {
//         Info << "error: " << sCNames_[sC] << endl;
        stateControllers_[sC]->calculateProperties();
    }
}


//- output -- call this function at the end of the MD time-step
void agentControllers::outputStateResults() 
{
    const Time& runTime = time_;

    if(runTime.outputTime())
    {
        // -- creating a set of directories in the current time directory
        {
            List<fileName> timePathNames(sCFixedPathNames_.size()); 
    
            if(nStateControllers_ > 0)
            {
                if(Pstream::master())
                {
                    // directory: case/<timeDir>/uniform
                    fileName uniformTimePath(runTime.path()/runTime.timeName()/"uniform");
                
                    if (!isDir(uniformTimePath))
                    {
                        mkDir(uniformTimePath);
                    }
        
        
                    if(stateControllers_.size() > 0)
                    {
                        // directory: case/<timeDir>/uniform/controllers
                        fileName controllersTimePath(uniformTimePath/"agentControllers");
    
                        if (!isDir(controllersTimePath))
                        {
                            mkDir(controllersTimePath);
                        }
    
    
                        forAll(stateControllers_, sC)
                        {
                            if(stateControllers_[sC]->writeInTimeDir())
                            {
                                // directory: case/<timeDir>/uniform/controllers/poly/<stateControllerModel>
                                fileName sCTimePath(controllersTimePath/sCNames_[sC]);
    
                                if(!isDir(sCTimePath))
                                {
                                    mkDir(sCTimePath);
                                }
    
                                //- creating directory for different zones but of the same model
                                const word& regionName = stateControllers_[sC]->regionName();
    
                                // directory: case/<timeDir>/uniform/controllers/poly/<stateControllerModel>/<cellZoneName>
                                fileName zoneTimePath(sCTimePath/regionName);
    
                                if (!isDir(zoneTimePath))
                                {
                                    mkDir(zoneTimePath);
                                }
    
                                timePathNames[sC] = zoneTimePath;
                            }
                        }
                    }
                }
            }

            // -- write out data (do not comment this out)
            forAll(stateControllers_, sC)
            {
                stateControllers_[sC]->output(sCFixedPathNames_[sC], timePathNames[sC]);
            }
        }
    }
}

const label& agentControllers::nStateControllers() const
{
    return nStateControllers_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
