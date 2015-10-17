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

#include "polyStateController.H"
#include "IFstream.H"
#include "graph.H"
#include "polyMoleculeCloud.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyStateController, 0);

defineRunTimeSelectionTable(polyStateController, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyStateController::polyStateController
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(molCloud.mesh())),
    molCloud_(molCloud),
    controllerDict_(dict.subDict("controllerProperties")),
	time_(t), 
    regionName_(controllerDict_.lookup("zoneName")),
    regionId_(-1),
    control_(true),
    controlInterForces_(false),
    readStateFromFile_(true),
    singleValueController_(false),
    density_(0.0),
    velocity_(vector::zero),
    temperature_(0.0),
    pressure_(0.0),
    strainRate_(tensor::zero),
    tempGradient_(vector::zero),
    fieldController_(false),
    densities_(),
    velocities_(),
    temperatures_(),
    pressures_(),
    writeInTimeDir_(true),
    writeInCase_(true)
{
    const cellZoneMesh& cellZones = mesh_.cellZones();
    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("polyStateController::polyStateController()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"controllersDict"
            << exit(FatalError);
    }

    control_ = Switch(controllerDict_.lookup("controlSwitch"));
    readStateFromFile_ = Switch(controllerDict_.lookup("readStateFromFile"));
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<polyStateController> polyStateController::New
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
{
    word polyStateControllerName
    (
        dict.lookup("stateControllerModel")
    );

    Info<< "Selecting polyStateController "
         << polyStateControllerName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyStateControllerName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "polyStateController::New(const dictionary&) : " << endl
            << "    unknown polyStateController type "
            << polyStateControllerName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<polyStateController>
	(
		cstrIter()(t, molCloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyStateController::~polyStateController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polyStateController::updateStateControllerProperties
(
    const dictionary& newDict
)
{
    controllerDict_ = newDict.subDict("controllerProperties");

    //- you can reset the controlling zone from here. This essentially
    //  means that the coupling zone can infact move arbitrarily. To make
    //  this happen we probably need to devise a technique for automatically
    //  changing the cellZone else where, and then calling this function to
    //  reset the controlling zone in which the controller operates in.

    if (controllerDict_.found("controlSwitch"))
    {
        control_ = Switch(controllerDict_.lookup("controlSwitch"));
    }

    if (controllerDict_.found("readStateFromFile"))
    {
        readStateFromFile_ = Switch(controllerDict_.lookup("readStateFromFile"));
    }
}

const labelList& polyStateController::controlZone() const
{
    return mesh_.cellZones()[regionId_];
}

const word& polyStateController::regionName() const
{
    return regionName_;
}

const bool& polyStateController::controlInterForces() const
{
    return controlInterForces_;
}

bool& polyStateController::controlInterForces()
{
    return controlInterForces_;
}

const scalar& polyStateController::density() const
{
    return density_;
}

scalar& polyStateController::density()
{
    return density_;
}

const vector& polyStateController::velocity() const
{
    return velocity_;
}

vector& polyStateController::velocity()
{
    return velocity_;
}

const scalar& polyStateController::temperature() const
{
    return temperature_;
}

scalar& polyStateController::temperature()
{
    return temperature_;
}

const scalar& polyStateController::pressure() const
{
    return pressure_;
}

scalar& polyStateController::pressure()
{
    return pressure_;
}

const tensor& polyStateController::strainRate() const
{
    return strainRate_;
}

tensor& polyStateController::strainRate()
{
    return strainRate_;
}

const vector& polyStateController::tempGradient() const
{
    return tempGradient_;
}

vector& polyStateController::tempGradient()
{
    return tempGradient_;
}

const scalarField& polyStateController::densityField() const
{
    return densities_;
}

scalarField& polyStateController::densityField()
{
    return densities_;
}

const vectorField& polyStateController::velocityField() const
{
    return velocities_;
}
vectorField& polyStateController::velocityField()
{
    return velocities_;
}

const scalarField& polyStateController::temperatureField() const
{
    return temperatures_;
}

scalarField& polyStateController::temperatureField()
{
    return temperatures_;
}


const scalarField& polyStateController::pressureField() const
{
    return pressures_;
}

scalarField& polyStateController::pressureField()
{
    return pressures_;
}


const bool& polyStateController::singleValueController() const
{
    return singleValueController_;
}

bool& polyStateController::singleValueController()
{
    return singleValueController_;
}

const bool& polyStateController::fieldController() const
{
    return fieldController_;
}

bool& polyStateController::fieldController()
{
    return fieldController_;
}

const bool& polyStateController::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& polyStateController::writeInCase() const
{
    return writeInCase_;
}



scalar polyStateController::avReqDensity()
{
    scalar totalDensity = 0.0;

    if(singleValueController_) 
    {
        totalDensity = density_;
    }
    else if(fieldController_)
    {
        label controlCells = controlZone().size();
    
        forAll(densities_, c)
        {
            totalDensity += densities_[c];
        }
    
        if (Pstream::parRun())
        {
            reduce(totalDensity, sumOp<scalar>());
            reduce(controlCells, sumOp<label>());
        }
    
        if(controlCells > 0)
        {
            totalDensity /= scalar(controlCells);
        }
    }

    return totalDensity;
}

vector polyStateController::avReqVelocity()
{
    vector totalVel = vector::zero;

    if(singleValueController_) 
    {
        totalVel = velocity_;
    }
    else if(fieldController_)
    {
        label controlCells = controlZone().size();
    
        forAll(velocities_, c)
        {
            totalVel += velocities_[c];
        }
    
        if (Pstream::parRun())
        {
            reduce(totalVel, sumOp<vector>());
            reduce(controlCells, sumOp<label>());
        }
    
        if(controlCells > 0)
        {
            totalVel /= scalar(controlCells);
        }
    }

    return totalVel;
}


scalar polyStateController::avReqTemperature()
{
    scalar totalTemp = 0.0;

    if(singleValueController_) 
    {
        totalTemp = temperature_;
    }
    else if(fieldController_)
    {
        label controlCells = controlZone().size();
    
        forAll(temperatures_, c)
        {
            totalTemp += temperatures_[c];
        }
    
        if (Pstream::parRun())
        {
            reduce(totalTemp, sumOp<scalar>());
            reduce(controlCells, sumOp<label>());
        }
    
        if(controlCells > 0)
        {
            totalTemp /= scalar(controlCells);
        }
    }
    return totalTemp;
}

scalar polyStateController::avReqPressure()
{
    scalar totalPressure = 0.0;

    if(singleValueController_) 
    {
        totalPressure = pressure_;
    }
    else if(fieldController_)
    {
        label controlCells = controlZone().size();
    
        forAll(pressures_, c)
        {
            totalPressure += pressures_[c];
        }
    
        if (Pstream::parRun())
        {
            reduce(totalPressure, sumOp<scalar>());
            reduce(controlCells, sumOp<label>());
        }
    
        if(controlCells > 0)
        {
            totalPressure /= scalar(controlCells);
        }
    }

    return totalPressure;
}

} // End namespace Foam

// ************************************************************************* //
