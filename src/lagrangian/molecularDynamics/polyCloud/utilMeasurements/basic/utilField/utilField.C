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

#include "utilField.H"



namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(utilField, 0);

defineRunTimeSelectionTable(utilField, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
utilField::utilField
(
    Time& t,
    const polyMesh& mesh,
    const reducedUnits& rU, 
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    time_(t),
    rU_(rU)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<utilField> utilField::New
(
    Time& t,
    const polyMesh& mesh,
    const reducedUnits& rU,
    const dictionary& dict
)
{
    word utilFieldName
    (
        dict.lookup("model")
    );

    Info<< "Selecting field: "
         << utilFieldName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(utilFieldName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "utilField::New(const dictionary&) : " << endl
            << "    unknown utilField type "
            << utilFieldName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<utilField>
	(
		cstrIter()(t, mesh, rU, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

utilField::~utilField()
{}


const fileName& utilField::outputPath() const
{
    return outputPath_;
}

fileName& utilField::outputPath()
{
    return outputPath_;
}

const fileName& utilField::inputPath() const
{
    return inputPath_;
}

fileName& utilField::inputPath()
{
    return inputPath_;
}


} // End namespace Foam

// ************************************************************************* //
