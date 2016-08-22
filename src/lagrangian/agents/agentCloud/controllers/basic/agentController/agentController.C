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

#include "agentController.H"
#include "IFstream.H"
#include "graph.H"
#include "agentCloud.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(agentController, 0);

defineRunTimeSelectionTable(agentController, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
agentController::agentController
(
    Time& t,
    agentCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(cloud.mesh())),
    cloud_(cloud),
//     controllerDict_(dict.subDict("controllerProperties")),
	time_(t), 
    regionName_(dict.lookup("zoneName")),
    regionId_(-1),
    control_(true),
    controlInterForces_(false),
    writeInTimeDir_(true),
    writeInCase_(true)
{
    const cellZoneMesh& cellZones = mesh_.cellZones();
    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("agentController::agentController()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"controllersDict"
            << exit(FatalError);
    }
    
    if(dict.found("control"))
    {    
        control_ = Switch(dict.lookup("control"));
    }
//     readStateFromFile_ = Switch(dict.lookup("readStateFromFile"));
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<agentController> agentController::New
(
    Time& t,
    agentCloud& cloud,
    const dictionary& dict
)
{
    word agentControllerName
    (
        dict.lookup("model")
    );

    Info<< "Selecting agentController "
         << agentControllerName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(agentControllerName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "agentController::New(const dictionary&) : " << endl
            << "    unknown agentController type "
            << agentControllerName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<agentController>
	(
		cstrIter()(t, cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentController::~agentController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// void agentController::updateStateControllerProperties
// (
//     const dictionary& newDict
// )
// {
// //     controllerDict_ = newDict.subDict("controllerProperties");
// 
//     //- you can reset the controlling zone from here. This essentially
//     //  means that the coupling zone can infact move arbitrarily. To make
//     //  this happen we probably need to devise a technique for automatically
//     //  changing the cellZone else where, and then calling this function to
//     //  reset the controlling zone in which the controller operates in.
// 
//     if (newDict.found("control"))
//     {
//         control_ = Switch(newDict.lookup("control"));
//     }
// /*
//     if (newDict.found("readStateFromFile"))
//     {
//         readStateFromFile_ = Switch(newDict.lookup("readStateFromFile"));
//     }*/
// }

const labelList& agentController::controlZone() const
{
    return mesh_.cellZones()[regionId_];
}

const word& agentController::regionName() const
{
    return regionName_;
}

const bool& agentController::controlInterForces() const
{
    return controlInterForces_;
}

bool& agentController::controlInterForces()
{
    return controlInterForces_;
}


const bool& agentController::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& agentController::writeInCase() const
{
    return writeInCase_;
}



} // End namespace Foam

// ************************************************************************* //
