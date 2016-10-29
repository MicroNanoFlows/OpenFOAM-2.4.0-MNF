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

#include "axialSelfDiffusion.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(axialSelfDiffusion, 0);

addToRunTimeSelectionTable(polyField, axialSelfDiffusion, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //






// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
axialSelfDiffusion::axialSelfDiffusion
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
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    useBoundBox_(false),
    unitVector_(propsDict_.lookup("unitVector")),    
    nSteps_(readLabel(propsDict_.lookup("nSteps")))  

{
    const cellZoneMesh& cellZones = mesh_.cellZones();

    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("axialSelfDiffusion::axialSelfDiffusion()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

   // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();
    
    nS_ = 0;
    nBatch_ = 0.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

axialSelfDiffusion::~axialSelfDiffusion()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void axialSelfDiffusion::createField()
{
    nMols_ = molCloud_.moleculeTracking().getMaxTrackingNumber();
    
    acf_.setSize(nMols_, 0.0);
    
    velocities_.clear();
    velocities_.setSize(nMols_);
    
    forAll(velocities_, i)
    {
        velocities_.setSize(nSteps_, 0.0);
        acf_.setSize(nSteps_, 0.0);
    }
    
    mols_.clear();
    mols_.setSize(nMols_, 0);

    // set initial velocities 
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    label actualNmols_ = 0;
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            label tN = molI->trackingNumber();                
            velocities_[tN][nS_] = molI->v() & unitVector_;
            mols_[tN] = 1;                
            actualNmols_ ++;
        }
    }    
    
    if(Pstream::parRun())
    {
        forAll(mols_, i)
        {
            reduce(mols_[i], sumOp<label>());
        }
        
        reduce(actualNmols_, sumOp<label>());
    }
    
    nS_++;
}


void axialSelfDiffusion::setVelocities()
{
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            label tN = molI->trackingNumber();                
            velocities_[tN][nS_] = molI->v() & unitVector_;
        }
    }     
}

void axialSelfDiffusion::calculateField()
{
    if(nS_ > nSteps_)
    {
        Info << "axialSelfDiffusion averaging " << endl;
        
        // parallel processing
        
        if(Pstream::parRun())
        {
            forAll(velocities_, i)
            {
                forAll(velocities_[i], j)
                {
                    reduce(velocities_[i][j], sumOp<scalar>());
                }
            }
        } 
    
        // calculate molecule based ACF 

        for (label i=0; i<nMols_; i++)
        {
            if(mols_[i] == 1)
            {
                forAll(velocities_[i], j)
                {
                    scalar u1=velocities_[i][j];
                    
                    forAll(velocities_[i], k)
                    {
                        scalar u1=velocities_[i][j];
                        scalar u11 = u1*u1;
                    }
                }
            }
        }
    
        // accumulate acf_
        nBatch_ += 1.0;        
        
        

        
        scalar integral = getIntegral();
        
        D_ = (1/(actualNmols_))*integral;
        
        nS_ = 0;

        
    }
    else
    {
        setVelocities();
        nS_++;
    }    
    
    
/*
    if(Pstream::parRun())
    {
        forAll(vDotV, i)
        {
            reduce(vDotV[i], sumOp<scalar>());
        }
    }*/
    
//     // sum vDotV over all molecules in the system
//     scalar sum = 0.0;
//     
//     forAll(vDotV, i)
//     {
//         sum += vDotV[i]*mols_[i];
//     }
    
    acfTime_[nS_] += sum;
   

}

// scalar axialSelfDiffusion::getIntegral()
// {
//     scalar timeIntegration = 0.0;
// 
//     const scalar& dt = time_.deltaT().value();
//     
//     if(((f.size() -1) % 2) == 0)// simpsons 1/3 rule
//     {
//         timeIntegration += f[0];
//         timeIntegration += f[f.size()-1];
// 
//         for (label i=1; i<f.size()-1; i++)
//         {
//             if((i % 2) == 0) // -even
//             {
//                 timeIntegration += 2.0*f[i];
//             }
//             else // odd
//             {
//                 timeIntegration += 4.0*f[i];
//             }
//         }
//         
//         timeIntegration *= dt/3.0;
// 
//     }
//     else // trapezoid rule
//     {
//         timeIntegration += f[0];
//         timeIntegration += f[f.size()-1];
// 
//         for (label i=1; i<f.size()-1; i++)
//         {
//             timeIntegration += 2.0*f[i];
//         }
// 
//         timeIntegration *= 0.5*dt;
//     }
//         
//     return timeIntegration;
// }

void axialSelfDiffusion::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {

        if(Pstream::master())
        {
            scalarField timeFieldI (1, runTime.timeOutputValue());
            scalarField D (1, D_);

            writeTimeData
            (
                casePath_,
                "axialSelfDiffusionCoeff_"+regionName_+"_"+fieldName_+"_D.xy",
                timeFieldI,
                D,
                true
            );
            
            scalarField timeFieldII (nSteps_, 0.0);
            scalarField acf (nSteps_, 0.0);
            
            const scalar& deltaT = time_.time().deltaT().value();
            
            forAll(timeFieldII, i)
            {
                timeFieldII[i]= deltaT*i;
                
                if(nAveragingSteps_ > 0)
                {
                    acf[i]=acfTime_[i]/nAveragingSteps_;
                }
            }
            
            writeTimeData
            (
                timePath_,
                "axialSelfDiffusionCoeff_"+regionName_+"_"+fieldName_+"_acf.xy",
                timeFieldII,
                acf
            );
        }
    }
}

void axialSelfDiffusion::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

void axialSelfDiffusion::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}


const propertyField& axialSelfDiffusion::fields() const
{
    return fields_;
}



} // End namespace Foam

// ************************************************************************* //
