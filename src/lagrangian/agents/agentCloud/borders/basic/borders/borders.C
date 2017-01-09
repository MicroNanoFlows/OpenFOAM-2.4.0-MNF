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

#include "borders.H"
#include "agentCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

borders::borders
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
    bordersDict_
    (
        IOobject
        (
            "bordersDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    
    borderList_(bordersDict_.lookup("borders")),
    borderNames_(borderList_.size()),
    borderIds_(borderList_.size()),
    borders_(borderList_.size())

{
    Info << nl << "Creating borders" << nl << endl;

    //- build borders

    if(borders_.size() > 0 )
    {
        forAll(borders_, i)
        {
            const entry& borderI = borderList_[i];
            const dictionary& bordersIDict = borderI.dict();
            
            borders_[i] = autoPtr<borderModel>
            (
                borderModel::New(time_, cloud_, bordersIDict)
            );
    
            borderNames_[i] = borders_[i]->type();
            borderIds_[i] = i;
        }
    }
    
    //test
    
    // making directory
/*    pathName_ = mesh_.time().path()/"borders";

    if(isDir(pathName_))
    {
        rmDir(pathName_);
    }     
    
    mkDir(pathName_);*/     
}

borders::~borders()
{}


void borders::initialConfig()
{
    forAll(borders_, i)
    {
//         borders_[i]->initialiseBorders();
        borders_[i]->initialConfiguration();  
    }
}

void borders::afterMove()
{
    forAll(borders_, i)
    {
        borders_[i]->afterMove();
    }
}

void borders::afterForce()
{
    forAll(borders_, i)
    {
        borders_[i]->afterForce();
    }
}

void borders::write()
{
    forAll(borders_, i)
    {
        borders_[i]->initialConfiguration();  
    }
}



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
