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

#include "scheduleModel.H"
#include "IFstream.H"
#include "graph.H"
#include "agentCloud.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(scheduleModel, 0);

defineRunTimeSelectionTable(scheduleModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
scheduleModel::scheduleModel
(
    Time& time,
    agentCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(cloud.mesh())),
    time_(time),
    cloud_(cloud),
    rU_(cloud_.redUnits())
{

}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<scheduleModel> scheduleModel::New
(
    Time& time,
    agentCloud& cloud,
    const dictionary& dict
)
{
    word scheduleModelName
    (
        dict.lookup("model")
    );

    Info<< "Selecting scheduleModel "
         << scheduleModelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(scheduleModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "scheduleModel::New(const dictionary&) : " << endl
            << "    unknown scheduleModel type "
            << scheduleModelName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<scheduleModel>
	(
		cstrIter()(time, cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

scheduleModel::~scheduleModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



} // End namespace Foam

// ************************************************************************* //
