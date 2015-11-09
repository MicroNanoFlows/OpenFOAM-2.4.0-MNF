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
#include "mathematicalConstants.H"

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
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    molIds_()
{

    writeInTimeDir_ = true;
    writeInCase_ = true;

//     singleValueController() = true;

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();
    
    // read in tracking numbers
//     trackingNumbers_ = List<label>(propsDict_.lookup("trackingNumbers"));
    trackingNumber_ = readLabel(propsDict_.lookup("trackingNumber"));

    centre_ = propsDict_.lookup("centre");
    a_ = readScalar(propsDict_.lookup("a"));
    b_ = readScalar(propsDict_.lookup("b"));
    deltaT_ = readScalar(propsDict_.lookup("deltaT"));    
    t_ = 0.0;
    
//     deltaT_ = time_.deltaT().value(); 
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tail::~tail()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tail::initialConfiguration()
{
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(mol().trackingNumber() == trackingNumber_)
        {
            Info << "centre y before = " << centre_.y() << endl;
            
            centre_.y() = mol().position().y()-b_;
            
            Info << "centre y after = " << centre_.y() << endl;
            
        }
    }    
    
}


void tail::controlBeforeVelocityI()
{
    Info << "tail: control" << endl;
    
    t_ += deltaT_;
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
//         if(findIndex(molIds_, mol().id()) != -1)
        {
//             Info << "mol().trackingNumber() = " << mol().trackingNumber() << endl;
            
            if(mol().trackingNumber() == trackingNumber_)
            {
                scalar xNew = -a_*cos(t_+(constant::mathematical::pi/2)) + centre_.x();
                scalar yNew = b_*sin(t_+(constant::mathematical::pi/2)) + centre_.y();
                
                mol().v().z() = 0;
                
                vector rNew = vector(xNew, yNew, 0.0);
                
                vector deltaR = rNew - mol().position();
                
                mol().a() = vector::zero;
                mol().v() = deltaR/deltaT_;
                Info << "position x = " << mol().position().x()
                     << ", position y = " << mol().position().y()
                     << endl;
            }
        }
    }
}

void tail::controlBeforeMove()
{}


void tail::controlBeforeForces()
{}

void tail::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void tail::controlAfterForces()
{}

void tail::controlAfterVelocityII()
{}

void tail::calculateProperties()
{}

void tail::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
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
