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

#include "dsmcCouplingController.H"

#include "IFstream.H"
#include "graph.H"
#include "dsmcCloud.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dsmcCouplingController, 0);

defineRunTimeSelectionTable(dsmcCouplingController, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcCouplingController::dsmcCouplingController
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict,
    couplingInterface1d& oneDInterfaces,
    couplingInterface2d& twoDInterfaces,
    couplingInterface3d& threeDInterfaces
)
:
    mesh_(cloud.mesh()),
    cloud_(cloud),
    rndGen_(cloud_.rndGen()),
    time_(t),
    regionNames_(dict.lookup("zoneNames")),
    regionIds_(),
    control_(true),
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

autoPtr<dsmcCouplingController> dsmcCouplingController::New
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict,
    couplingInterface1d& oneDInterfaces,
    couplingInterface2d& twoDInterfaces,
    couplingInterface3d& threeDInterfaces
)
{
    word dsmcCouplingControllerName
    (
        dict.lookup("couplingControllerModel")
    );

    Info<< "Selecting dsmcCouplingController "
         << dsmcCouplingControllerName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(dsmcCouplingControllerName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "dsmcCouplingController::New(const dictionary&) : " << endl
            << "    unknown dsmcCouplingController type "
            << dsmcCouplingControllerName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<dsmcCouplingController>
	(
		cstrIter()(t, cloud, dict, oneDInterfaces, twoDInterfaces, threeDInterfaces)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcCouplingController::~dsmcCouplingController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcCouplingController::updateCouplingControllerProperties
(
    const dictionary& newDict
)
{
    //- you can reset the controlling zone from here. This essentially
    //  means that the coupling zone can in fact move arbitrarily. To make
    //  this happen we probably need to devise a technique for automatically
    //  changing the cellZone else where, and then calling this function to
    //  reset the controlling zone in which the controller operates in.

    if (newDict.found("control"))
    {
        control_ = Switch(newDict.lookup("control"));
    }
}

const labelList& dsmcCouplingController::controlZone(label regionID) const
{
    return mesh_.cellZones()[regionID];
}

const List<word>& dsmcCouplingController::regionNames() const
{
    return regionNames_;
}

const List<label>& dsmcCouplingController::regionIds() const
{
    return regionIds_;
}

const bool& dsmcCouplingController::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& dsmcCouplingController::writeInCase() const
{
    return writeInCase_;
}

} // End namespace Foam

// ************************************************************************* //
