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

#include "polyCouplingPatch.H"

#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyCouplingPatch, 0);

addToRunTimeSelectionTable(polyPatchBoundary, polyCouplingPatch, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyCouplingPatch::polyCouplingPatch
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyPatchBoundary(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    elapsedTime_(0.0),
    writeInterval_(readScalar(t.controlDict().lookup("writeInterval"))),
    startTime_(t.startTime().value()),
    molFlux_(0.0),
    massFlux_(0.0),
    cumulMolFlux_(0.0),
    cumulMassFlux_(0.0)
{
    writeInTimeDir_ = false;
    writeInCase_ = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyCouplingPatch::~polyCouplingPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void polyCouplingPatch::initialConfiguration()
{}

void polyCouplingPatch::calculateProperties()
{}

void polyCouplingPatch::controlMol
(
    polyMolecule& mol,
    polyMolecule::trackingData& td
)
{}

void polyCouplingPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}

void polyCouplingPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
}

} // End namespace Foam

// ************************************************************************* //
