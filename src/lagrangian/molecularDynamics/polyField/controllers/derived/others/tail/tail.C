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

#include "tail.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(tail, 0);

addToRunTimeSelectionTable(polyStateController, tail, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
tail::tail
(
    Time& t,
//     const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, /*mesh,*/ molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    molIds_()
{

    writeInTimeDir_ = true;
    writeInCase_ = true;

    singleValueController() = true;

    molIds_.clear();

    selectIds ids
    (
        molCloud_.pot(),
        propsDict_
    );

    molIds_ = ids.molIds();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tail::~tail()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tail::initialConfiguration()
{}

void tail::calculateProperties()
{}

void tail::controlMolsBeg()
{}

void tail::controlBeforeForces()
{}

void tail::controlMols()
{
    Info << "tail: control" << endl;

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            const polyMolecule::constantProperties& constProp = molCloud_.constProps(mol().id());
            const scalar& massI = constProp.mass();

            mol().a() += force/massI;

            force_ += force;
        }
    }
}

void tail::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

// void tail::controlDuringForces
// (
//     polyMolecule* molReal,
//     polyReferredMolecule* molRef
// )
// {}

void tail::controlMolsEnd()
{}

void tail::output
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
//                 "tail_force.xyz",
//                 timeField,
//                 forces,
//                 true
//             );            
//         }
    }
}

void tail::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");


}

} // End namespace Foam

// ************************************************************************* //
