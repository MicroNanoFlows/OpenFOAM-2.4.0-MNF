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

#include "dsmcCouplingPatch.H"

#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcCouplingPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcCouplingPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcCouplingPatch::dsmcCouplingPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    allSpecies_(false),
    typeIds_()
{
    measurePropertiesAtWall_ = false;
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcCouplingPatch::~dsmcCouplingPatch()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcCouplingPatch::initialConfiguration()
{}

void dsmcCouplingPatch::calculateProperties()
{}

void dsmcCouplingPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
{
    td.keepParticle = false;

    //Add to list of coupled parcels before it is deleted
    dsmcParcel *parcelPtr = &p;
    cloud_.insertCoupledParcel(parcelPtr);
}

void dsmcCouplingPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}

void dsmcCouplingPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    setProperties();
}

void dsmcCouplingPatch::setProperties()
{
    if (propsDict_.found("allSpecies"))
    {
        allSpecies_ = Switch(propsDict_.lookup("allSpecies"));
    }

    if(!allSpecies_)
    {
        const List<word> molecules (propsDict_.lookup("typeIds"));

        DynamicList<word> moleculesReduced(0);

        forAll(molecules, i)
        {
            const word& moleculeName(molecules[i]);

            if(findIndex(moleculesReduced, moleculeName) == -1)
            {
                moleculesReduced.append(moleculeName);
            }
        }

        moleculesReduced.shrink();

        typeIds_.clear();

        typeIds_.setSize(moleculesReduced.size(), -1);

        forAll(moleculesReduced, i)
        {
            const word& moleculeName(moleculesReduced[i]);

            label typeId(findIndex(cloud_.typeIdList(), moleculeName));

            if(typeId == -1)
            {
                FatalErrorIn("dsmcCouplingPatch::dsmcCouplingPatch()")
                    << "Cannot find typeId: " << moleculeName << nl << "in: "
                    << mesh_.time().system()/"boundariesDict"
                    << exit(FatalError);
            }

            typeIds_[i] = typeId;
        }
    }
}


} // End namespace Foam

// ************************************************************************* //
