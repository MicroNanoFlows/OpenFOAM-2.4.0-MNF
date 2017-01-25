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

#include "agentZoneProperties.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agentZoneProperties, 0);
addToRunTimeSelectionTable(agentMeasurement, agentZoneProperties, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


void agentZoneProperties::setBoundBoxes()
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


// Construct from components
agentZoneProperties::agentZoneProperties
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
    boxes_(),
    totalArea_(0.0),
    N_(0.0),
    mass_(0.0),
    mom_(vector::zero),
    vel_(vector::zero),
    kE_(0.0),    
    timeIndex_(0),
    NField_(1, 0.0),
    rhoNField_(1, 0.0),
    rhoMField_(1, 0.0),
    UField_(1, vector::zero),
    momField_(1, vector::zero),
    kEField_(1, 0.0),
    nTimeSteps_(0.0),
    resetAtOutput_(true)  
{
   
    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();

    setBoundBoxes();
    
    forAll(boxes_, b)
    {
        totalArea_ += (boxes_[b].span() & vector(1, 0 ,0) ) * (boxes_[b].span() & vector(0, 1 ,0) );
    }
    
    Info << "total area = " <<  totalArea_ << endl;
    
    resetAtOutput_ = Switch(propsDict_.lookup("resetAtOutput"));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentZoneProperties::~agentZoneProperties()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void agentZoneProperties::createField()
{
}

void agentZoneProperties::calculateField()
{
    nTimeSteps_ += 1.0;
    
    {
        scalar N = 0.0;
        scalar mass = 0.0;
        vector mom = vector::zero;
        vector vel = vector::zero;
        scalar kE = 0.0;  
    
        IDLList<agent>::iterator mol(cloud_.begin());
        
        for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
        {
            if(findIndex(agentIds_, mol().id()) != -1)
            {
                forAll(boxes_, b)
                {
                    if(boxes_[b].contains(mol().position()))
                    {
                        N += 1.0;
                        mass += mol().mass();
                        mom += mol().v()*mol().mass();
                        vel += mol().v();
                        kE += 0.5*magSqr(mol().v());
                    }
                }
            }
        }
            
        N_ += N;
        mass_ += mass;
        mom_ += mom;
        vel_ += vel;
        kE_ += kE;        
    }
    
    if(time_.outputTime())
    {
        scalar N = N_;
        scalar mass = mass_;
        vector mom = mom_;
        vector vel = vel_;
        scalar kE = kE_;         
            
        if(Pstream::parRun())
        {
            reduce(N, sumOp<scalar>());
            reduce(mass, sumOp<scalar>());
            reduce(mom, sumOp<vector>());
            reduce(vel, sumOp<vector>());
            reduce(kE, sumOp<scalar>());
        }         
        
        
        
        NField_[timeIndex_] = N/nTimeSteps_;
        rhoNField_[timeIndex_] = N/(totalArea_*nTimeSteps_);
        rhoMField_[timeIndex_] = mass/(totalArea_*nTimeSteps_);
        UField_[timeIndex_] = vel/nTimeSteps_;
        kEField_[timeIndex_] = kE/nTimeSteps_;
        momField_[timeIndex_] = mom/nTimeSteps_;
        
        if(resetAtOutput_)
        {
            nTimeSteps_ = 0.0;
            N_ = 0.0;
            mass_ = 0.0;
            mom_ = vector::zero;
            vel_ = vector::zero;
            kE_ = 0.0;
        }
    }
}

void agentZoneProperties::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            scalarField timeField(1, runTime.timeOutputValue());
            
            writeTimeData
            (
                casePath_,
                "zone_bb_"+fieldName_+"_N.xy",
                timeField,
                NField_,
                true
            );
            
            writeTimeData
            (
                casePath_,
                "zone_bb_"+fieldName_+"_rhoN.xy",
                timeField,
                rhoNField_,
                true
            );

            writeTimeData
            (
                casePath_,
                "zone_bb_"+fieldName_+"_rhoM.xy",
                timeField,
                rhoMField_,
                true
            );

            writeTimeData
            (
                casePath_,
                "zone_bb_"+fieldName_+"_U.xyz",
                timeField,
                UField_,
                true
            );

            writeTimeData
            (
                casePath_,
                "zone_bb_"+fieldName_+"_mom.xyz",
                timeField,
                momField_,
                true
            );            
            
            writeTimeData
            (
                casePath_,
                "zone_bb_"+fieldName_+"_kE.xy",
                timeField,
                kEField_,
                true
            );


        }
    }    
}



void agentZoneProperties::measureDuringForceComputation
(
    agent* molI,
    agent* molJ
)
{}

void agentZoneProperties::measureDuringForceComputationSite
(
    agent* molI,
    agent* molJ,
    label sI,
    label sJ
)
{}

// const propertyField& agentZoneProperties::fields() const
// {
//     return fields_;
// }

} // End namespace Foam

// ************************************************************************* //
