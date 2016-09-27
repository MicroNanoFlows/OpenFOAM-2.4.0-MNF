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

#include "schedules.H"
#include "agentCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

schedules::schedules
(
    Time& t,
    const polyMesh& mesh,
    agentCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    time_(t),
    cloud_(cloud),
    cP_(cloud_.cP()),
    schedulesDict_
    (
        IOobject
        (
            "schedulesDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    
    scheduleList_(schedulesDict_.lookup("schedules")),
    scheduleNames_(scheduleList_.size()),
    scheduleIds_(scheduleList_.size()),
    schedules_(scheduleList_.size())

{
    Info << nl << "Creating schedules" << nl << endl;

    //- build schedules

    if(schedules_.size() > 0 )
    {
        forAll(schedules_, i)
        {
            const entry& scheduleI = scheduleList_[i];
            const dictionary& schedulesIDict = scheduleI.dict();
            
            schedules_[i] = autoPtr<scheduleModel>
            (
                scheduleModel::New(time_, cloud_, schedulesIDict)
            );
    
            scheduleNames_[i] = schedules_[i]->type();
            scheduleIds_[i] = i;
        }
    }
    
    //test
    
    // making directory
/*    pathName_ = mesh_.time().path()/"schedules";

    if(isDir(pathName_))
    {
        rmDir(pathName_);
    }     
    
    mkDir(pathName_);*/     
}

schedules::~schedules()
{}

/*
void schedules::initialConfig()
{
    forAll(schedules_, i)
    {
        schedules_[i]->initialiseBorders();
        schedules_[i]->initialConfiguration();  
    }
}

void schedules::afterMove()
{
    forAll(schedules_, i)
    {
        schedules_[i]->afterMove();
    }
}

void schedules::afterForce()
{
    forAll(schedules_, i)
    {
        schedules_[i]->afterForce();
    }
}

void schedules::write()
{
    
}*/



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
