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
    agentWallForce

Description

\*----------------------------------------------------------------------------*/

#include "agentWallForce.H"
#include "graph.H"
#include "IFstream.H"
#include "OFstream.H"
#include "agentCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(agentWallForce, 0);

defineRunTimeSelectionTable(agentWallForce, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

agentWallForce::agentWallForce
(
    Time& time,
    const dictionary& dict
)
:
    time_(time),
    spaceVarying_(false),
    timeVarying_(false)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<agentWallForce> agentWallForce::New
(
    Time& time,
    const dictionary& dict
)
{
    word agentWallForceName
    (
        dict.lookup("model")
    );

    Info<< "Selecting wall-force model "
         << agentWallForceName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(agentWallForceName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "agentWallForce::New(const dictionary&) : " << endl
            << "    unknown agentWallForce type "
            << agentWallForceName
            << ", constructor not in hash table" << endl << endl
            << "    Valid agentWallForce types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<agentWallForce>
    (
        cstrIter()(time, dict)
    );
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentWallForce::~agentWallForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool agentWallForce::spaceVarying()
{
    return spaceVarying_;
}

bool agentWallForce::timeVarying()
{
    return timeVarying_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
