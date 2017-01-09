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

#include "agentMeasurement.H"



namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(agentMeasurement, 0);

defineRunTimeSelectionTable(agentMeasurement, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
agentMeasurement::agentMeasurement
(
    Time& t,
    const polyMesh& mesh,
    agentCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    time_(t),
    casePath_(),
    timePath_(),
    measureInterForces_(false),
    measureInterForcesSites_(false)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<agentMeasurement> agentMeasurement::New
(
    Time& t,
    const polyMesh& mesh,
    agentCloud& cloud,
    const dictionary& dict
)
{
    word agentMeasurementName
    (
        dict.lookup("model")
    );

    Info<< "Selecting field: "
         << agentMeasurementName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(agentMeasurementName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "agentMeasurement::New(const dictionary&) : " << endl
            << "    unknown agentMeasurement type "
            << agentMeasurementName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<agentMeasurement>
	(
		cstrIter()(t, mesh, cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentMeasurement::~agentMeasurement()
{}


const fileName& agentMeasurement::casePath() const
{
    return casePath_;
}

fileName& agentMeasurement::casePath()
{
    return casePath_;
}

const fileName& agentMeasurement::timePath() const
{
    return timePath_;
}

fileName& agentMeasurement::timePath()
{
    return timePath_;
}

const bool& agentMeasurement::measureInterForces() const
{
    return measureInterForces_;
}

bool& agentMeasurement::measureInterForces()
{
    return measureInterForces_;
}

const bool& agentMeasurement::measureInterForcesSites() const
{
    return measureInterForcesSites_;
}

bool& agentMeasurement::measureInterForcesSites()
{
    return measureInterForcesSites_;
}
} // End namespace Foam

// ************************************************************************* //
