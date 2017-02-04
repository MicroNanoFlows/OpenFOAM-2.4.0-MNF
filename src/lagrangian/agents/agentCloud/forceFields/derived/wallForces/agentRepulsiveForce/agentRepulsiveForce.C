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
    agentRepulsiveForce

Description

\*----------------------------------------------------------------------------*/

#include "agentRepulsiveForce.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(agentRepulsiveForce, 0);

addToRunTimeSelectionTable(agentWallForce, agentRepulsiveForce, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
agentRepulsiveForce::agentRepulsiveForce
(
    Time& time,
    const dictionary& dict
)
:
    agentWallForce(time, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    force_(propsDict_.lookup("force"))
{
    timeVarying_ = true;
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentRepulsiveForce::~agentRepulsiveForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector agentRepulsiveForce::force(const vector& position)
{
    return vector::zero;
}

void agentRepulsiveForce::updateForce()
{

}

vector agentRepulsiveForce::force(const scalar& time)
{
    return force_;
}

vector agentRepulsiveForce::force(const scalar& time, const vector& position)
{
    return force_;
}

void agentRepulsiveForce::write
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
