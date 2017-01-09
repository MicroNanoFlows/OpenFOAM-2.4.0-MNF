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

#include "basicOutputs.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(basicOutputs, 0);

addToRunTimeSelectionTable(agentMeasurement, basicOutputs, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
basicOutputs::basicOutputs
(
    Time& t,
    const polyMesh& mesh,
    agentCloud& molCloud,
    const dictionary& dict
)
:
    agentMeasurement(t, mesh, molCloud, dict),
//     fields_(t, mesh, "dummy"),
    accumulatedTotalLinearMomentum_(vector::zero),
    accumulatedTotalMass_(0.0),
    accumulatedTotalLinearKE_(0.0),
    accumulatedTotalPE_(0.0),
    accumulatedNMols_(0.0),

    nAvTimeSteps_(0.0)
{
    const scalarField& cellVols = mesh_.cellVolumes();
    
    meshVolume_ = sum(cellVols);
    
    if (Pstream::parRun())
    {
        reduce(meshVolume_, sumOp<scalar>());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

basicOutputs::~basicOutputs()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void basicOutputs::createField()
{}


void basicOutputs::calculateField()
{
    nAvTimeSteps_ += 1.0;
    
    vector singleStepTotalLinearMomentum(vector::zero);
    scalar singleStepMaxVelocityMag = 0.0;
    scalar singleStepTotalMass = 0.0;
    scalar singleStepTotalLinearKE = 0.0;
    scalar singleStepTotalPE = 0.0;

    label singleStepNMols = cloud_.size();

    {
        IDLList<agent>::iterator mol(cloud_.begin());
    
        for
        (
            mol = cloud_.begin();
            mol != cloud_.end();
            ++mol
        )
        {
//             label molId = mol().id();
            scalar molMass = mol().mass();
            singleStepTotalMass += molMass;
        }
    
        for
        (
            mol = cloud_.begin();
            mol != cloud_.end();
            ++mol
        )
        {
//             label molId = mol().id();
    
            scalar molMass = mol().mass();

            const vector& molV(mol().v());

            singleStepTotalLinearMomentum += molV * molMass;

            if(mag(molV) > singleStepMaxVelocityMag)
            {
                singleStepMaxVelocityMag = mag(molV);
            }
    
            singleStepTotalLinearKE += 0.5*molMass*magSqr(molV);
            singleStepTotalPE += mol().potentialEnergy();
        }
    }

    if (Pstream::parRun())
    {
        reduce(singleStepTotalLinearMomentum, sumOp<vector>());
        reduce(singleStepMaxVelocityMag, maxOp<scalar>());
        reduce(singleStepTotalMass, sumOp<scalar>());
        reduce(singleStepTotalLinearKE, sumOp<scalar>());
        reduce(singleStepTotalPE, sumOp<scalar>());
        reduce(singleStepNMols, sumOp<label>());
    }

    if (singleStepNMols)
    {
        Info<< "Number of molecules in system = "
            << singleStepNMols << nl
            << "Overall number density = "
            << singleStepNMols/meshVolume_ << nl
            << "Overall mass density = "
            << singleStepTotalMass/meshVolume_ << nl
            << "Average linear momentum per molecule = "
            << singleStepTotalLinearMomentum/singleStepNMols << ' '
            << mag(singleStepTotalLinearMomentum)/singleStepNMols << nl
            << "Maximum |velocity| = "
            << singleStepMaxVelocityMag << nl
            << "Average linear KE per molecule = "
            << singleStepTotalLinearKE/singleStepNMols << nl
            << "Average PE per molecule = "
            << singleStepTotalPE/singleStepNMols << nl
            << "Average TE per molecule = "
            <<
            (
                singleStepTotalLinearKE
            + singleStepTotalPE
            )
            /singleStepNMols
            << endl;
    }
    else
    {
        Info<< "No molecules in system" << endl;
    }

    accumulatedTotalLinearMomentum_ += singleStepTotalLinearMomentum;
    accumulatedTotalMass_ += singleStepTotalMass;
    accumulatedTotalLinearKE_ += singleStepTotalLinearKE;
    accumulatedTotalPE_ += singleStepTotalPE;
    accumulatedNMols_ += singleStepNMols;

}

void basicOutputs::writeField()
{
    const Time& runTime = time_.time();

    if (runTime.outputTime())
    {
//     	const scalar& nAvTimeSteps = nAvTimeSteps_; 
// 
//         if (accumulatedNMols_)
//         {
// 
//         }
//         else
//         {
//             Info<< "Not averaging temperature and pressure: "
//                 << "no molecules in system" << endl;
//         }

        //-reset
        accumulatedTotalLinearMomentum_ = vector::zero;
        accumulatedTotalMass_ = 0.0;
        accumulatedTotalLinearKE_ = 0.0;
        accumulatedTotalPE_ = 0.0;
        accumulatedNMols_ = 0.0;
        nAvTimeSteps_ = 0.0;
    }
}

void basicOutputs::measureDuringForceComputation
(
    agent* molI,
    agent* molJ
)
{}

void basicOutputs::measureDuringForceComputationSite
(
    agent* molI,
    agent* molJ,
    label sI,
    label sJ
)
{}

// const propertyField& basicOutputs::fields() const
// {
//     return  fields_;
// }

} // End namespace Foam

// ************************************************************************* //
