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

#include "polyBoundedForce.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyBoundedForce, 0);

addToRunTimeSelectionTable(polyStateController, polyBoundedForce, dictionary);


void polyBoundedForce::setBoundBoxes()
{
 
    PtrList<entry> boxList(propsDict_.lookup("boxes"));

    boxes_.setSize(boxList.size());

    forAll(boxList, b)
    {
        const entry& boxI = boxList[b];
        const dictionary& dict = boxI.dict();

        vector startPoint = dict.lookup("startPoint");
        vector endPoint = dict.lookup("endPoint");
        boxes_[b].resetBoundedBox(startPoint, endPoint);
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyBoundedForce::polyBoundedForce
(
    Time& t,
//     const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, /*mesh,*/ molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    model_(),
    molIds_(),
    nTimeSteps_(0.0),
    force_(vector::zero),
    timeVarying_(false)
    
{

    writeInTimeDir_ = true;
    writeInCase_ = true;

    singleValueController() = true;

    model_ = autoPtr<gravityForce>
    (
        gravityForce::New(t, propsDict_)
    );

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    readProperties();

    setBoundBoxes();
    
    // incorporate the ability for this model to handle both
    // time-varying and spatially varying flows
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyBoundedForce::~polyBoundedForce()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyBoundedForce::initialConfiguration()
{}

void polyBoundedForce::calculateProperties()
{}

void polyBoundedForce::controlMolsBeg()
{}

void polyBoundedForce::controlBeforeForces()
{}

void polyBoundedForce::controlMols()
{
    // time changes in forces
    model_->updateForce();

    // - if control switch is on
    if(control_)
    {
        Info << "polyBoundedForce: control" << endl;

        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                forAll(boxes_, b)
                {
                    if(boxes_[b].contains(mol().position()))
                    {
                        const scalar& massI = molCloud_.cP().mass(mol().id());
                        
                        vector force = model_->force(mol().position());

                        mol().a() += force/massI;

                        force_ += force;
                    }
                }
            }
        }
        
        nTimeSteps_ += 1.0;
    }
}

void polyBoundedForce::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

// void polyBoundedForce::controlDuringForces
// (
//     polyMolecule* molReal,
//     polyReferredMolecule* molRef
// )
// {}

void polyBoundedForce::controlMolsEnd()
{}

void polyBoundedForce::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    model_->write(fixedPathName, timePath);

    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
//         vector force = force_;
// 
//         if(Pstream::parRun())
//         {
//             reduce(force, sumOp<vector>());
//         }
// 
//         if(Pstream::master())
//         {
//             vectorField forces(1, force/nTimeSteps_);
//             
//             scalarField timeField(1, time_.time().timeOutputValue());
//    
//             writeTimeData
//             (
//                 fixedPathName,
//                 "polyBoundedForce_force.xyz",
//                 timeField,
//                 forces,
//                 true
//             );            
//         }
    }
}

void polyBoundedForce::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    model_->updateProperties(propsDict_);

    readProperties();
}


void polyBoundedForce::readProperties()
{

}



} // End namespace Foam

// ************************************************************************* //