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

#include "routes.H"
#include "agentCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

routes::routes
(
    Time& t,
    const polyMesh& mesh,
    agentCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    time_(t),
    cloud_(cloud),
    routesDict_
    (
        IOobject
        (
            "routesDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    
    routeList_(routesDict_.lookup("routes")),
    routeNames_(routeList_.size()),
    routeIds_(routeList_.size()),
    routes_(routeList_.size())

{
    Info << nl << "Creating routes" << nl << endl;

    //- build routes

    if(routes_.size() > 0 )
    {
        forAll(routes_, i)
        {
            const entry& routeI = routeList_[i];
            const dictionary& routesIDict = routeI.dict();
            
            routes_[i] = autoPtr<routeModel>
            (
                scheduleModel::New(time_, cloud_, routesIDict)
            );
    
            routeNames_[i] = routes_[i]->type();
            routeIds_[i] = i;
        }
    }
}

routes::~routes()
{}


void routes::initialConfig()
{
    forAll(routes_, i)
    {
        routes_[i]->initialConfiguration();  
    }
}


/*
void routes::afterMove()
{
    forAll(routes_, i)
    {
        routes_[i]->afterMove();
    }
}

void routes::afterForce()
{
    forAll(routes_, i)
    {
        routes_[i]->afterForce();
    }
}

void routes::write()
{
    
}*/



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
