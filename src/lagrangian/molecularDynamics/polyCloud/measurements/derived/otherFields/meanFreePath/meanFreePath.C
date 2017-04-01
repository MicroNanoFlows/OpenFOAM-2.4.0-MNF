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

#include "meanFreePath.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(meanFreePath, 0);

addToRunTimeSelectionTable(polyField, meanFreePath, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //






// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
meanFreePath::meanFreePath
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
    molIds_()
{
    dCol_ = readScalar(propsDict_.lookup("collisionDistance"));
    
    
   // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();
    
    
    deltaT_ = time_.deltaT().value();    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meanFreePath::~meanFreePath()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void meanFreePath::createField()
{
    label nTNs = molCloud_.moleculeTracking().getMaxTrackingNumber();
    
    Info << "number of tracking numbers = " << nTNs << endl;
    listSize_ = nTNs;
    
    freePaths_.setSize(listSize_, 0.0);
    Rold_.setSize(listSize_, 0.0);
    collectorOfFreePaths_.setSize(listSize_);
    collided_.setSize(listSize_, false);
    setRolds();
}

void meanFreePath::setRolds()
{
    Rold_ = 0.0;
    
    {
        IDLList<agent>::iterator mol(molCloud_.begin());
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            label tNI = mol().trackingNumber();
            
            Rold_[tNI] = mol().R();   
        }
    }    
    
    // parallel processing
    
    forAll(Rold_, i)
    {
        if(Pstream::parRun())
        {
            reduce(Rold_[i], sumOp<scalar>());
        }
    }
}

void meanFreePath::calculateField()
{
    // free paths 
    
    List<scalar> freePaths.setSize(listSize_, 0.0);
    
    {
        IDLList<agent>::iterator mol(molCloud_.begin());
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            label tNI = mol().trackingNumber();
            
            if(mol().R() > dCol_) // free path
            {
                freePaths[tNI] += mol().v()*deltaT_;
            }
        }
    }
    
    // parallel comms 

    forAll(freePaths, i)
    {
        if(Pstream::parRun())
        {
            reduce(freePaths[i], sumOp<scalar>());
        }
    }
    
    freePaths_ += freePaths;
    
    // check for collisions
    collided_ = false;
    
    {
        IDLList<agent>::iterator mol(molCloud_.begin());
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            label tNI = mol().trackingNumber();
            
            if( (mol().R() < dCol_) && (Rold_[tNI] > dCol_) ) // collision
            {
                collided_[tNI] = true;
            }
        }
    } 
    
    setRolds();
    
    
    Info << "av number of collisions per molecule = " << << endl;
}

void meanFreePath::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {

        if(Pstream::master())
        {

            massFlux_.shrink();
            scalarField timeField (massFlux_.size(), 0.0);
            scalarField massFlux (massFlux_.size(), 0.0);
            
            massFlux.transfer(massFlux_);
            massFlux_.clear();
            
            const scalar& deltaT = time_.time().deltaT().value();
            
            forAll(timeField, i)
            {
                timeField[timeField.size()-i-1]=runTime.timeOutputValue()-(deltaT*i);
            }
            
            writeTimeData
            (
                casePath_,
                "fluxZone_"+regionName_+"_"+fieldName_+"_M.xy",
                timeField,
                massFlux,
                true
            );


            const reducedUnits& rU = molCloud_.redUnits();
    
            if(rU.outputSIUnits())
            {
                writeTimeData
                (
                    casePath_,
                    "fluxZone_"+regionName_+"_"+fieldName_+"_M_SI.xy",
                    timeField*rU.refTime(),
                    massFlux*rU.refMassFlux(),
                    true
                );
            }
        }
    }
}

void meanFreePath::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void meanFreePath::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}


const propertyField& meanFreePath::fields() const
{
    return fields_;
}



} // End namespace Foam

// ************************************************************************* //
