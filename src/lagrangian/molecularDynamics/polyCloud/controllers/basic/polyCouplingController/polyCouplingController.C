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

#include "polyCouplingController.H"

#include "IFstream.H"
#include "graph.H"
#include "polyMoleculeCloud.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyCouplingController, 0);
defineRunTimeSelectionTable(polyCouplingController, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyCouplingController::polyCouplingController
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict,
    couplingInterface2d& twoDInterfaces,
    couplingInterface3d& threeDInterfaces
)
:
    mesh_(refCast<const fvMesh>(molCloud.mesh())),
    molCloud_(molCloud),
    time_(t),
    regionNames_(dict.lookup("zoneNames")),
    regionIds_(),
    control_(true),
    controlInterForces_(false),
    writeInTimeDir_(true),
    writeInCase_(true)
{
    const cellZoneMesh& cellZones = mesh_.cellZones();

    forAll(regionNames_, regions)
    {
      label lclRegionId = cellZones.findZoneID(regionNames_[regions]);

      if(lclRegionId == -1)
      {
          FatalErrorIn("dsmcCouplingController::dsmcCouplingController()")
              << "Cannot find region: " << regionNames_[regions] << nl << "in: "
              << time_.time().system()/"controllersDict"
              << exit(FatalError);
      }

      regionIds_.append(lclRegionId);
    }

    if(dict.found("control"))
    {
        control_ = Switch(dict.lookup("control"));
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<polyCouplingController> polyCouplingController::New
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict,
    couplingInterface2d& twoDInterfaces,
    couplingInterface3d& threeDInterfaces
)
{
    word polyCouplingControllerName
    (
        dict.lookup("couplingControllerModel")
    );

    Info<< "Selecting polyCouplingController "
         << polyCouplingControllerName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyCouplingControllerName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "polyCouplingController::New(const dictionary&) : " << endl
            << "    unknown polyCouplingController type "
            << polyCouplingControllerName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<polyCouplingController>
    (
        cstrIter()(t, molCloud, dict, twoDInterfaces, threeDInterfaces)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyCouplingController::~polyCouplingController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polyCouplingController::updateCouplingControllerProperties
(
    const dictionary& newDict
)
{
    //- you can reset the controlling zone from here. This essentially
    //  means that the coupling zone can infact move arbitrarily. To make
    //  this happen we probably need to devise a technique for automatically
    //  changing the cellZone else where, and then calling this function to
    //  reset the controlling zone in which the controller operates in.

    if (newDict.found("control"))
    {
        control_ = Switch(newDict.lookup("control"));
    }
}

const labelList& polyCouplingController::controlZone(label regionID) const
{
    return mesh_.cellZones()[regionID];
}

const List<word>& polyCouplingController::regionNames() const
{
    return regionNames_;
}

const List<label>& polyCouplingController::regionIds() const
{
    return regionIds_;
}

const bool& polyCouplingController::controlInterForces() const
{
    return controlInterForces_;
}

bool& polyCouplingController::controlInterForces()
{
    return controlInterForces_;
}

const bool& polyCouplingController::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& polyCouplingController::writeInCase() const
{
    return writeInCase_;
}

} // End namespace Foam

// ************************************************************************* //
