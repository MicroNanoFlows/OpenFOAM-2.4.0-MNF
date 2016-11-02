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
//     regionName_(propsDict_.lookup("zoneName")),
//     regionId_(-1),
    useBoundBox_(false),
    unitVector_(propsDict_.lookup("unitVector")),    
    nSteps_(readLabel(propsDict_.lookup("nSteps")))  

{
//     const cellZoneMesh& cellZones = mesh_.cellZones();
// 
//     regionId_ = cellZones.findZoneID(regionName_);
// 
//     if(regionId_ == -1)
//     {
//         FatalErrorIn("axialSelfDiffusion::axialSelfDiffusion()")
//             << "Cannot find region: " << regionName_ << nl << "in: "
//             << time_.time().system()/"fieldPropertiesDict"
//             << exit(FatalError);
//     }

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
    
    tNaddress_.setSize(nMols_);
    
    
    label actualNmols_ = 0;
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    label actualNmols_ = 0;
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            label tN = mol().trackingNumber();
            tNaddress_[tN] = actualNmols_;
            actualNmols_ ++;
        }
    }
    
    //- parallel communication
    if(Pstream::parRun())
    {
        //-sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);

                    toNeighbour << mols << mass << mom 
                                << kE << angularKeSum << dof << kineticTensor 
                                << virialTensor;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                scalarField molsProc;
                scalarField massProc;
                vectorField momProc;
                scalarField kEProc;
                scalarField angularKeSumProc;
                scalarField dofProc;
                tensorField kineticTensorProc;
                tensorField virialTensorProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour  >> molsProc >> massProc >> momProc 
                                >> kEProc >> angularKeSumProc >> dofProc >> kineticTensorProc
                                >> virialTensorProc;
                }
    
                mols += molsProc;
                mass += massProc;
                mom += momProc;
                kE += kEProc;
                angularKeSum += angularKeSumProc;
                dof += dofProc;
                kineticTensor += kineticTensorProc;
                virialTensor += virialTensorProc;
            }
        }
    }
    
    
    acf_.setSize(nMols_);
    
    velocities_.clear();
    velocities_.setSize(nMols_);
    
    forAll(velocities_, i)
    {
        velocities_[i].setSize(nSteps_, 0.0);
        acf_[i].setSize(nSteps_, 0.0);
    }
    
    mols_.clear();
    mols_.setSize(nMols_, 0);

    ACF_.setSize(nSteps_, 0.0);
    
    // set initial velocities 
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            label tN = mol().trackingNumber();                
            velocities_[tN][nS_] = mol().v() & unitVector_;
            mols_[tN] = 1;                
        }
    }    
    
    
    nS_++;
    
    
        FatalErrorIn("axialSelfDiffusion::axialSelfDiffusion()")
            << "Fatal error test "
            << exit(FatalError);    
}


void axialSelfDiffusion::setVelocities()
{
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            label tN = mol().trackingNumber();                
            velocities_[tN][nS_] = mol().v() & unitVector_;
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
        
        label T = nSteps_;
        
        
        for (label i=0; i<nMols_; i++)
        {
            if(mols_[i] == 1)
            {
                for (label k=0; k<T; k++)
                {
                    scalar uSum = 0.0;
                    
                    for (label t=0; t<T-k; t++)
                    {
                        uSum += velocities_[i][t]*velocities_[i][t+k];
                    }
                
                    acf_[i][k] += uSum/(T-1);
                }
            }
        }
    

        nBatch_ += 1.0;        
        
        // the acf summed across molecule batches and all molecules
        
        for (label i=0; i<nMols_; i++)
        {
            if(mols_[i] == 1)
            {
                for (label k=0; k<T; k++)
                {
                    ACF_[k] += acf_[i][k]/(nBatch_*scalar(actualNmols_));
                }
            }
        }
        
        
        // integrate ACF_ 
        
        nS_ = 0;
    }
    else
    {
        setVelocities();
        nS_++;
    }    
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
            scalarField timeFieldII (nSteps_, 0.0);
            scalarField acf (nSteps_, 0.0);
            
            const scalar& deltaT = time_.time().deltaT().value();
            
            forAll(timeFieldII, i)
            {
                timeFieldII[i]= deltaT*i;
                acf[i]=ACF_[i];
            }
            
            writeTimeData
            (
                timePath_,
                "ACF_"+fieldName_+".xy",
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
