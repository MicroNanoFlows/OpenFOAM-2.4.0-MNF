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


void axialSelfDiffusion::setBoundBox
(
    const dictionary& propsDict,
    boundedBox& bb,
    const word& name 
)
{
    const dictionary& dict(propsDict.subDict(name));
    
    vector startPoint = dict.lookup("startPoint");
    vector endPoint = dict.lookup("endPoint");
    bb.resetBoundedBox(startPoint, endPoint);
}




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


    //-set the total volume
//     const labelList& cells = cellZones[regionId_];
// 
//     forAll(cells, c)
//     {
//         const label& cellI = cells[c];
//         totalVolume_ += mesh_.cellVolumes()[cellI];
//     }
// 
//     if (Pstream::parRun())
//     {
//         reduce(totalVolume_, sumOp<scalar>());
//     }
    
    if (propsDict_.found("useBoundBox"))
    {
        useBoundBox_ = Switch(propsDict_.lookup("useBoundBox"));
        
        if(useBoundBox_)
        {
            setBoundBox(propsDict_, bb_, "samplingRegion");
        }
    }    
    
    nS_ = 0;
    nAveragingSteps_ = 0.0;
    acfTime_.setSize(nSteps_, 0.0);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

axialSelfDiffusion::~axialSelfDiffusion()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void axialSelfDiffusion::createField()
{
    nMols_ = molCloud_.moleculeTracking().getMaxTrackingNumber();
    
    
    initialVelocities_.clear();
    initialVelocities_.setSize(nMols_, 0.0);

    
    mols_.clear();
    mols_.setSize(nMols_, 0);

    // set initial velocities 
    
    const List< DynamicList<polyMolecule*> >& cellOccupancy
            = molCloud_.cellOccupancy();
    
    const labelList& cells = mesh_.cellZones()[regionId_];
    
    forAll(cells, c)
    {
        const label& cellI = cells[c];
        const List<polyMolecule*>& molsInCell = cellOccupancy[cellI];
        
        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];
            
            if(findIndex(molIds_, molI->id()) != -1)
            {
                label tN = molI->trackingNumber();
                
//                 const scalar& massI = molCloud_.cP().mass(molI->id());
                
                if(useBoundBox_)
                {
                    if(bb_.contains(molI->position()))
                    {
                        initialVelocities_[tN] = molI->v() & unitVector_;
                        mols_[tN] = 1;
                    }
                }
                else
                {
                    initialVelocities_[tN] = molI->v() & unitVector_;
                    mols_[tN] = 1;
                }
            }
        }
    }
    
    if(Pstream::parRun())
    {
        forAll(initialVelocities_, i)
        {
            reduce(initialVelocities_[i], sumOp<scalar>());
            reduce(mols_[i], sumOp<label>());
        }
    }    
        
}

void axialSelfDiffusion::calculateField()
{
    // set current velocities 
    
    const List< DynamicList<polyMolecule*> >& cellOccupancy
            = molCloud_.cellOccupancy();
    
    const labelList& cells = mesh_.cellZones()[regionId_];
    
    List<scalar> vDotV(nMols_, 0.0);
        
    forAll(cells, c)
    {
        const label& cellI = cells[c];
        const List<polyMolecule*>& molsInCell = cellOccupancy[cellI];
        
        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];
            
            if(findIndex(molIds_, molI->id()) != -1)
            {
                label tN = molI->trackingNumber();
                
//                 const scalar& massI = molCloud_.cP().mass(molI->id());
                
                if(useBoundBox_)
                {
                    if(bb_.contains(molI->position()))
                    {
                        vDotV[tN] = (molI->v() & unitVector_) * initialVelocities_[tN];
                    }
                }
                else
                {
                    vDotV[tN] =  (molI->v() & unitVector_) * initialVelocities_[tN];
                }
            }
        }
    }    

    if(Pstream::parRun())
    {
        forAll(vDotV, i)
        {
            reduce(vDotV[i], sumOp<scalar>());
        }
    }
    
    // sum vDotV over all molecules in the system
    scalar sum = 0.0;
    
    forAll(vDotV, i)
    {
        sum += vDotV[i]*mols_[i];
    }
    
    acfTime_[nS_] += sum;
   
    if(nS_ >= nSteps_)
    {
        Info << "axialSelfDiffusion averaging " << endl; 
        nS_ = 0;
        nAveragingSteps_ += 1.0;
        
        List<scalar> acfTime = acfTime_;
        
        forAll(acfTime, i)
        {
            acfTime[i] /= nAveragingSteps_;
        }
        
        scalar integral = getIntegral(acfTime);
        
        D_ = (1/(3.0*nMols_))*integral;
        
        //reset field
        createField(); 
    }
    else
    {
        nS_++;    
    }
}

scalar axialSelfDiffusion::getIntegral(const List<scalar>& f)
{
    scalar timeIntegration = 0.0;

    const scalar& dt = time_.deltaT().value();
    
    if(((f.size() -1) % 2) == 0)// simpsons 1/3 rule
    {
        timeIntegration += f[0];
        timeIntegration += f[f.size()-1];

        for (label i=1; i<f.size()-1; i++)
        {
            if((i % 2) == 0) // -even
            {
                timeIntegration += 2.0*f[i];
            }
            else // odd
            {
                timeIntegration += 4.0*f[i];
            }
        }
        
        timeIntegration *= dt/3.0;

    }
    else // trapezoid rule
    {
        timeIntegration += f[0];
        timeIntegration += f[f.size()-1];

        for (label i=1; i<f.size()-1; i++)
        {
            timeIntegration += 2.0*f[i];
        }

        timeIntegration *= 0.5*dt;
    }
        
    return timeIntegration;
}

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
