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

#include "agent1DBins.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agent1DBins, 0);
addToRunTimeSelectionTable(agentMeasurement, agent1DBins, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from components
agent1DBins::agent1DBins
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
//     totalArea_(0.0),
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
    
    // create bin model
    binModel_ = autoPtr<binModel>
    (
        binModel::New(mesh, propsDict_)
    );
    
    nBins_ = binModel_->nBins();
    
    N_.setSize(nBins_, 0.0);
    mass_.setSize(nBins_, 0.0);
    mom_.setSize(nBins_, vector::zero);
    vel_.setSize(nBins_, vector::zero);
    kE_.setSize(nBins_, 0.0);

    NField_.setSize(nBins_, 0.0);
    rhoNField_.setSize(nBins_, 0.0);
    rhoMField_.setSize(nBins_, 0.0);
    UField_.setSize(nBins_, vector::zero);
    momField_.setSize(nBins_, vector::zero);
    kEField_.setSize(nBins_, 0.0);

    resetAtOutput_ = Switch(propsDict_.lookup("resetAtOutput"));
    
    boundedBox bbMesh( mesh_.bounds().min(), mesh_.bounds().max() );
    
       
    dZ_ = bbMesh.span() & vector(0, 0, 1);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agent1DBins::~agent1DBins()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void agent1DBins::createField()
{
}

void agent1DBins::calculateField()
{
    nTimeSteps_ += 1.0;
    
    {
        scalarField N(nBins_, 0.0);
        scalarField mass(nBins_, 0.0);
        vectorField mom(nBins_, vector::zero);
        vectorField vel(nBins_, vector::zero);
        scalarField kE(nBins_, 0.0);  
    
        IDLList<agent>::iterator mol(cloud_.begin());
        
        for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
        {
            if(findIndex(agentIds_, mol().id()) != -1)
            {
                label n = binModel_->isPointWithinBin(mol().position(), mol().cell());

                if(n != -1)
                {         
                    N[n] += 1.0;
                    mass[n] += mol().mass();
                    mom[n] += mol().v()*mol().mass();
                    vel[n] += mol().v();
                    kE[n] += 0.5*magSqr(mol().v());
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
        scalarField N = N_;
        scalarField mass = mass_;
        vectorField mom = mom_;
        vectorField vel = vel_;
        scalarField kE = kE_;         
            
        if(Pstream::parRun())
        {
            forAll(N, n)
            {
                reduce(N[n], sumOp<scalar>());
                reduce(mass[n], sumOp<scalar>());
                reduce(mom[n], sumOp<vector>());
                reduce(vel[n], sumOp<vector>());
                reduce(kE[n], sumOp<scalar>());
            }
        }         

        forAll(N, n)
        {
            scalar area = binModel_->binVolume(n)*dZ_;

            NField_[n] = N[n]/nTimeSteps_;
            rhoNField_[n] = N[n]/(area*nTimeSteps_);
            rhoMField_[n] = mass[n]/(area*nTimeSteps_);
            UField_[n] = vel[n]/nTimeSteps_;
            kEField_[n] = kE[n]/nTimeSteps_;
            momField_[n] = mom[n]/nTimeSteps_;
        }
        
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

void agent1DBins::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            scalarField bins = binModel_->binPositions();
             
            writeTimeData
            (
                timePath_,
                "bins_1D_"+fieldName_+"_N.xy",
                bins,
                NField_
            );
            
            writeTimeData
            (
                timePath_,
                "bins_1D_"+fieldName_+"_rhoN.xy",
                bins,
                rhoNField_
            );
            
            writeTimeData
            (
                timePath_,
                "bins_1D_"+fieldName_+"_rhoM.xy",
                bins,
                rhoMField_
            );            
            
            writeTimeData
            (
                timePath_,
                "bins_1D_"+fieldName_+"_U.xyz",
                bins,
                UField_
            ); 
            
            writeTimeData
            (
                timePath_,
                "bins_1D_"+fieldName_+"_mom.xyz",
                bins,
                momField_
            );

            writeTimeData
            (
                timePath_,
                "bins_1D_"+fieldName_+"_kE.xyz",
                bins,
                kEField_
            );
        }
    }    
}



void agent1DBins::measureDuringForceComputation
(
    agent* molI,
    agent* molJ
)
{}

void agent1DBins::measureDuringForceComputationSite
(
    agent* molI,
    agent* molJ,
    label sI,
    label sJ
)
{}

// const propertyField& agent1DBins::fields() const
// {
//     return fields_;
// }

} // End namespace Foam

// ************************************************************************* //
