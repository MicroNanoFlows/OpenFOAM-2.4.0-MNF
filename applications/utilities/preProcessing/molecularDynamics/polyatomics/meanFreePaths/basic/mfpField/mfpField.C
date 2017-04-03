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

#include "mfpField.H"



namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mfpField, 0);

defineRunTimeSelectionTable(mfpField, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
mfpField::mfpField
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

autoPtr<mfpField> mfpField::New
(
    Time& t,
    const polyMesh& mesh,
    const reducedUnits& rU,
    const dictionary& dict
)
{
    word mfpFieldName
    (
        dict.lookup("model")
    );

    Info<< "Selecting field: "
         << mfpFieldName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(mfpFieldName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "mfpField::New(const dictionary&) : " << endl
            << "    unknown mfpField type "
            << mfpFieldName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<mfpField>
	(
		cstrIter()(t, mesh, rU, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mfpField::~mfpField()
{}


const fileName& mfpField::outputPath() const
{
    return outputPath_;
}

fileName& mfpField::outputPath()
{
    return outputPath_;
}

const fileName& mfpField::inputPath() const
{
    return inputPath_;
}

fileName& mfpField::inputPath()
{
    return inputPath_;
}


} // End namespace Foam

// ************************************************************************* //
