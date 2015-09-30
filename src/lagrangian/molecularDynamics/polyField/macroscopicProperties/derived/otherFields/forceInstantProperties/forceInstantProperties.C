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

#include "forceInstantProperties.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(forceInstantProperties, 0);

addToRunTimeSelectionTable(polyField, forceInstantProperties, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void forceInstantProperties::setBoundBoxes()
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
forceInstantProperties::forceInstantProperties
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fields_(t, mesh, "dummy"),
    fieldName_(propsDict_.lookup("fieldName")),
    boxes_(),
    molIds_()
{

        // build bound boxes
    setBoundBoxes();

    // choose molecule ids to sample
    molIds_.clear();

    selectIds ids
    (
        molCloud_.pot(),
        propsDict_
    );

    molIds_ = ids.molIds();

    forceField_.clear();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forceInstantProperties::~forceInstantProperties()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void forceInstantProperties::createField()
{
 
}

void forceInstantProperties::calculateField()
{
    vector force = vector::zero;

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            forAll(boxes_, b)
            {
                if(boxes_[b].contains(mol().position()))
                {
                    const polyMolecule::constantProperties& constProp = molCloud_.constProps(mol().id());
                    force += constProp.mass()*mol().a();
                }
            }
        }
    }    
    
    // - parallel processing
    if(Pstream::parRun())
    {
        reduce(force, sumOp<vector>());
    }
    
    forceField_.append(force);
}


void forceInstantProperties::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            forceField_.shrink();

            scalarField timeField (forceField_.size(), 0.0);
            vectorField force (forceField_.size(), vector::zero);
            
            force.transfer(forceField_);
            forceField_.clear();
            
            const scalar& deltaT = time_.time().deltaT().value();
            
            forAll(timeField, i)
            {
                timeField[timeField.size()-i-1]=runTime.timeOutputValue()-(deltaT*i);
            }
            
            writeTimeData
            (
                casePath_,
                "force_instant_"+fieldName_+".xyz",
                timeField,
                force,
                true
            );
        }
    }
}
void forceInstantProperties::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

void forceInstantProperties::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}

const propertyField& forceInstantProperties::fields() const
{
    return fields_;
}

// void forceInstantProperties::updateProperties(const dictionary& newDict)
// {
//     //- the main properties should be updated first
//     updateBasicFieldProperties(newDict);

// }

} // End namespace Foam

// ************************************************************************* //
