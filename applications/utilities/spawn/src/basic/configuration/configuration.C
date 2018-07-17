/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "configuration.H"
#include "IFstream.H"
#include "graph.H"


namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(configuration, 0);
defineRunTimeSelectionTable(configuration, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
configuration::configuration
(
    molecules& cloud,
    const dictionary& dict
)
:
    cloud_(cloud),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<configuration> configuration::New
(

    molecules& cloud,
    const dictionary& dict
)
{
    word configurationName
    (
        dict.lookup("model")
    );

    Info<< "Creating new configuration model: " << configurationName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(configurationName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "configuration::New(const dictionary&) : " << endl
            << "    unknown configuration type "
            << configurationName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<configuration>
    (
        cstrIter()(cloud, dict)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

configuration::~configuration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void configuration::insertAtom
(
    const vector& position,
    const label type,
    const scalar& charge
)
{
    cloud_.positions().append(position);
    cloud_.types().append(type);
    cloud_.charges().append(charge);
}

void configuration::deleteAtom
(
    label id
)
{
    
}
/*

void configuration::checkForOverlaps()
{
    
    DynamicList<vector> positionsNew;
    
    scalar tolerance = 0.1;
    
    forAll (positions, i)
    {
        const vector& rI = positions[i];
        
        bool overlapping = false;
        
        forAll (positionsNew, j)
        {
            const vector& rJ = positionsNew[j];
            scalar rMag = mag(rI - rJ);
        
            if (rMag < tolerance)
            {
                overlapping = true;
            }
        }
    
        if (!overlapping)
        {
            positionsNew.append(rI);
        }
    }
}*/

/*
label configuration::mols() const
{
    return nMolsAdded_;
}*/

} // End namespace Foam

// ************************************************************************* //
