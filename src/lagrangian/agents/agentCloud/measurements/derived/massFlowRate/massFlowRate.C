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

#include "massFlowRate.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(massFlowRate, 0);
addToRunTimeSelectionTable(agentMeasurement, massFlowRate, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
massFlowRate::massFlowRate
(
    Time& t,
    const polyMesh& mesh,
    agentCloud& cloud,
    const dictionary& dict
)
:
    agentMeasurement(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    agentIds_(),
    fieldName_(propsDict_.lookup("fieldName")),
    unitVector_(propsDict_.lookup("flowDirection"))
{
    unitVector_ /= mag(unitVector_);
    
    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();
    
    massFlux_.clear();
    
    boundedBox bbMesh( mesh_.bounds().min(), mesh_.bounds().max() );
    
    length_ = bbMesh.span() & unitVector_;
    
    Info << "length = " << length_ << endl;
    
    nTimeSteps_ = 0.0;
    
    massFlowRate_ = 0.0;
    
    instant_ = false;
    
    if (propsDict_.found("instant"))
    {
        instant_ = Switch(propsDict_.lookup("instant"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

massFlowRate::~massFlowRate()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void massFlowRate::createField()
{
}

void massFlowRate::calculateField()
{
    IDLList<agent>::iterator mol(cloud_.begin());
    
    vector mom = vector::zero;
    
    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {
        if(findIndex(agentIds_, mol().id()) != -1)
        {
            mom += mol().v()*mol().mass();
        }
    }
    
    if(Pstream::parRun())
    {
        reduce(mom, sumOp<vector>());
    }    

    scalar massFlux = (mom & unitVector_)/(length_);
    
    massFlux_.append(massFlux);
 
    massFlowRate_ += massFlux;
    
    nTimeSteps_ += 1.0;
}

void massFlowRate::writeField()
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
            
            // output instantaneous file (may be a large file)
            if(instant_)
            {
                writeTimeData
                (
                    casePath_,
                    "massflowRate_"+fieldName_+"_instant.xy",
                    timeField,
                    massFlux,
                    true
                );
            }

            scalarField timeFieldII (1, 0.0);
            scalarField massFlowAv (1, 0.0);            
            
            timeFieldII[0] = runTime.timeOutputValue();
            
            if(nTimeSteps_ > 0)
            {
                massFlowAv[0] = massFlowRate_/nTimeSteps_;
            }
            
            // cumulative averaging
            writeTimeData
            (
                casePath_,
                "massflowRate_"+fieldName_+"_average.xy",
                timeFieldII,
                massFlowAv,
                true
            );            
        }
    }    
}



void massFlowRate::measureDuringForceComputation
(
    agent* molI,
    agent* molJ
)
{}

void massFlowRate::measureDuringForceComputationSite
(
    agent* molI,
    agent* molJ,
    label sI,
    label sJ
)
{}

// const propertyField& massFlowRate::fields() const
// {
//     return fields_;
// }

} // End namespace Foam

// ************************************************************************* //
