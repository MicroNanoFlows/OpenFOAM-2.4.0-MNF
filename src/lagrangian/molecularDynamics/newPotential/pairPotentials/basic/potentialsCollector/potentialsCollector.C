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

#include "potentialsCollector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// constructor
potentialsCollector::potentialsCollector
(
    const polyMesh& mesh,
    const reducedUnits& redUnits
)
:
    mesh_(mesh),
    potentialsDict_
    (
        IOobject
        (
            "potentialsDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    potentialsList_(potentialsDict_.lookup("pairs")),
    names_(potentialsList_.size()),
    ids_(potentialsList_.size()),
    potentials_(potentialsList_.size())
{
    if(potentialsList_.size() > 0 )
    {
        Info << nl << "Creating potentials: " << nl << endl;
    
        forAll(potentialsList_, i)
        {
            const entry& potI = potentialsList_[i];
            const dictionary& potIDict = potI.dict();
    
            potentials_[i] = autoPtr<potentialModel>
            (
                potentialModel::New(mesh, redUnits, potIDict)
            );
    
            names_[i] = potentials_[i]->type();
            ids_[i] = i;
        }
        
//         testPairPotentials();        
    }
    else
    {
         FatalErrorIn("potentialsCollector.C") << nl
                << " Something went wrong. You have no potentials in your system/potentialsDict"
                << nl << nl << "Check that the format is as follows: " 
                << " pairs " << nl
                << " ( " << nl
                << "    // enter your pair potentials here
                << " ); " 
                << nl << abort(FatalError);       
    }
}


potentialsCollector::~potentialsCollector()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// void potentialsCollector::testPairPotentials()
// {
//     Info << nl << "Initialising the measurement fields" << nl << endl;
// 
//     forAll(potentials_, f)
//     {
//         potentials_[f]->header();
//     }
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
