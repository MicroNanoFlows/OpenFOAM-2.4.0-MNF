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

#include "polyDelFromMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyDelFromMesh, 0);

addToRunTimeSelectionTable(polyMolsToDeleteModel, polyDelFromMesh, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyDelFromMesh::polyDelFromMesh
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyMolsToDeleteModel(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    oneSpecie_(false),
    molId_(-1)
{
    if (propsDict_.found("oneSpecie"))
    {
        oneSpecie_ = Switch(propsDict_.lookup("oneSpecie"));

        if (oneSpecie_)
        {
            const List<word>& idList(molCloud_.pot().idList());
            const word molId = propsDict_.lookup("molId");
            molId_ = findIndex(idList, molId);

            if(molId_ == -1)
            {
                FatalErrorIn("polyDelFromMesh::polyDelFromMesh()")
                    << "Cannot find molId: " << molId << nl << "in: "
                    << mesh_.time().system()/"molsToDeleteDict"
                    << exit(FatalError);
            }
        }
    }

    findMolsToDel();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyDelFromMesh::~polyDelFromMesh()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void polyDelFromMesh::findMolsToDel()
{
    DynamicList<polyMolecule*> molsToDel;

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

//     label counter = 0;
    label initialSize = molCloud_.size();

    for
    (
        mol = molCloud_.begin();
        mol != molCloud_.end();
        ++mol
    )
    {
        label molId = mol().id();

        if(!oneSpecie_ || (molId == molId_))
        {
            polyMolecule* molI = &mol();
            molsToDel.append(molI);
        }
    }

    //molsToDel.shrink();

    forAll(molsToDel, m)
    {
        deleteMolFromMoleculeCloud(*molsToDel[m]);
    }

    label molsKept = initialSize - molsToDel.size();

    Info<< tab << " initial polyMolecules: " <<  initialSize 
        << ", polyMolecules kept: " <<  molsKept
        << ", polyMolecules removed: " << molsToDel.size() 
        << endl;


    // as a precaution: rebuild cell occupancy
    molCloud_.rebuildCellOccupancy();
}


} // End namespace Foam

// ************************************************************************* //
