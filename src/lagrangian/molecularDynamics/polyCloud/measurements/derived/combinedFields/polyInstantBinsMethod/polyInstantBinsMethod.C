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

#include "polyInstantBinsMethod.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyInstantBinsMethod, 0);
addToRunTimeSelectionTable(polyField, polyInstantBinsMethod, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyInstantBinsMethod::polyInstantBinsMethod
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
    binModel_(),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    fieldName_(propsDict_.lookup("fieldName")),

    molIds_(),

    mols_(),
    mass_(),
    mom_(),
    velocityB_(),
    kE_(),
    angularKeSum_(),
    dof_(),
    kineticTensor_(),
    virialTensor_(),

    molsV_(),
    massV_(),
    momV_(),
    velocity_(),
    angularSpeed_(),
    angularVelocity_(),

    N_(),
    rhoN_(),
    rhoM_(),
    USAM_(),
    UCAM_(),
    T_(),
    stress_(),
    p_(),
    outputField_(4, true)
//     nAvTimeSteps_(0.0),
//     resetAtOutput_(true)
{
 
    
    const cellZoneMesh& cellZones = mesh_.cellZones();
    
    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("polyInstantBinsMethod::polyInstantBinsMethod()")
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

    
    // create bin model
    binModel_ = autoPtr<binModel>
    (
        binModel::New(mesh, propsDict_)
    );

    const label& nBins = binModel_->nBins();

    mols_.setSize(nBins, 0.0);
    mass_.setSize(nBins, 0.0);
    mom_.setSize(nBins, vector::zero);
    velocityB_.setSize(nBins, vector::zero);
    kE_.setSize(nBins, 0.0);
    angularKeSum_.setSize(nBins, 0.0);
    dof_.setSize(nBins, 0.0);
    kineticTensor_.setSize(nBins, tensor::zero);
    virialTensor_.setSize(nBins, tensor::zero);

    molsV_.setSize(nBins, 0.0);
    massV_.setSize(nBins, 0.0);
    momV_.setSize(nBins, vector::zero);
    velocity_.setSize(nBins, vector::zero);
    angularSpeed_.setSize(nBins, vector::zero);
    angularVelocity_.setSize(nBins, vector::zero);

    N_.setSize(nBins);
    rhoN_.setSize(nBins);
    rhoM_.setSize(nBins);
    USAM_.setSize(nBins);
    UCAM_.setSize(nBins);
    T_.setSize(nBins);
    stress_.setSize(nBins);
    p_.setSize(nBins);

    nBins_ = nBins;
    
    // choice of measurement property to output

    if (propsDict_.found("outputProperties"))
    {

        const List<word> measurements (propsDict_.lookup("outputProperties"));

        DynamicList<word> propertyNames(0);

        forAll(measurements, i)
        {
            const word& propertyName(measurements[i]);
    
            if(findIndex(propertyNames, propertyName) == -1)
            {
                propertyNames.append(propertyName);
            }
        }

        //propertyNames.shrink();

        if(findIndex(propertyNames, "density") == -1)
        {
            outputField_[0] = false;
        }

        if(findIndex(propertyNames, "velocity") == -1)
        {
            outputField_[1] = false;
        }

        if(findIndex(propertyNames, "temperature") == -1)
        {
            outputField_[2] = false;
        }

        if(findIndex(propertyNames, "pressure") == -1)
        {
            outputField_[3] = false;
        }

        // check for user:
        forAll(propertyNames, i)
        {
            const word& propertyName(propertyNames[i]);

            if
            (
                (propertyName != "density") &&
                (propertyName != "velocity") &&
                (propertyName != "temperature") &&
                (propertyName != "pressure")
            )
            {    
                FatalErrorIn("polyInstantBinsMethod::polyInstantBinsMethod()")
                    << "Cannot find measurement property: " << propertyName
                    << nl << "in: "
                    << time_.time().system()/"fieldPropertiesDict"
                    << exit(FatalError);            
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyInstantBinsMethod::~polyInstantBinsMethod()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyInstantBinsMethod::createField()
{
    Info << "Initialising polyInstantBinsMethod fields" << endl;

    const List< DynamicList<polyMolecule*> >& cellOccupancy
        = molCloud_.cellOccupancy();

    const labelList& cells = mesh_.cellZones()[regionId_];

    scalarField mols(mass_.size(), 0.0);
    scalarField mass(mass_.size(), 0.0);
    vectorField mom(mom_.size(), vector::zero);
    vectorField angularSpeed(mom_.size(), vector::zero);

    forAll(cells, c)
    {
        const label& cellI = cells[c];
        const List<polyMolecule*>& molsInCell = cellOccupancy[cellI];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            const vector& rI = molI->position();

            label n = binModel_->isPointWithinBin(rI, cellI);

            if(n != -1)
            {
                if(findIndex(molIds_, molI->id()) != -1)
                {
//                     const polyMolecule::constantProperties& constProp 
//                                 = molCloud_.constProps(molI->id());

                    mols[n] += 1.0;

                    const scalar& massI = molCloud_.cP().mass(molI->id());
                    
                    mass[n] += massI;

                    mom[n] += massI*molI->v();

                    const diagTensor& molMoI(molCloud_.cP().momentOfInertia(molI->id()));

                    // angular speed 
                    const vector& molOmega(inv(molMoI) & molI->pi());

                    angularSpeed[n] += molOmega;
                }
            }
        }
    }


    // parallel processing
    if(Pstream::parRun())
    {
        // sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << mols << mass << mom << angularSpeed;
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
                vectorField angularSpeedProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> molsProc >> massProc >> momProc >> angularSpeedProc;
                }

                mols += molsProc;
                mass += massProc;
                mom += momProc;
                angularSpeed += angularSpeedProc;
            }
        }
    }

    forAll(velocity_, n)
    {
        if(mass[n] > 0.0)
        {
            angularVelocity_[n] = angularSpeed[n]/mols[n];
            velocity_[n] = mom[n]/mass[n];
        }
    }
}


void polyInstantBinsMethod::calculateField()
{
    scalarField mols(nBins_, 0.0);
    scalarField mass(nBins_, 0.0);
    vectorField mom(nBins_, vector::zero);
    vectorField angularSpeed(nBins_, vector::zero);

    
    forAll(mesh_.cellZones()[regionId_], c)
    {
        const label& cellI = mesh_.cellZones()[regionId_][c];
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            const vector& rI = molI->position();

            label n = binModel_->isPointWithinBin(rI, cellI);

            if(n != -1)
            {
                if(findIndex(molIds_, molI->id()) != -1)
                {
                    mols[n] += 1.0;
                    
                    const scalar& massI = molCloud_.cP().mass(molI->id());
                    
                    mass[n] += massI;

                    mom[n] += massI*molI->v();

                    const diagTensor& molMoI(molCloud_.cP().momentOfInertia(molI->id()));

                    // angular speed 
                    const vector& molOmega(inv(molMoI) & molI->pi());

                    angularSpeed[n] += molOmega;
                }
            }
        }
    }

    //- parallel communication
    if(Pstream::parRun())
    {
        forAll(mass, i)
        {
            forAll(mass[i], j)
            {
                reduce(mols[i][j], sumOp<scalar>());
                reduce(mass[i][j], sumOp<scalar>());
                reduce(mom[i][j], sumOp<vector>());
                reduce(angularSpeed[i][j], sumOp<vector>());
            }
        }        
    }
    
    vectorField velocity(nBins_, vector::zero);
    vectorField angularVelocity(nBins_, vector::zero);
    
    forAll(velocity, n)
    {
        if(mass[n] > 0.0)
        {
            velocity[n] = mom[n]/mass[n];
            angularVelocity[n] = angularSpeed[n]/mols[n];
        }
    }
    
    scalarField dof(nBins_, 0.0);
    scalarField kE(nBins_, 0.0);
    scalarField angularKeSum(nBins_, 0.0);
    tensorField kineticTensor(nBins_, tensor::zero);
    tensorField virialTensor(nBins_, tensor::zero);
   
    forAll(mesh_.cellZones()[regionId_], c) 
    {
        const label& cellI = mesh_.cellZones()[regionId_][c];
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            const vector& rI = molI->position();

            label n = binModel_->isPointWithinBin(rI, cellI);

            if(n != -1)
            {
                if(findIndex(molIds_, molI->id()) != -1)
                {
                    const scalar& massI = molCloud_.cP().mass(molI->id());

                    kE[n] += 0.5*massI*magSqr(molI->v() - velocity_[n]);

                    dof[n] += molCloud_.cP().degreesOfFreedom(molI->id());

                    const diagTensor& molMoI(molCloud_.cP().momentOfInertia(molI->id()));

                    // angular speed 
                    const vector& molOmega(inv(molMoI) & molI->pi());

                    angularKeSum[n] += 0.5*(molOmega & molMoI & molOmega);

                    kineticTensor[n] += ( massI*(molI->v() - velocity_[n])*(molI->v() - velocity_[n]) ) 
                                            + ( ((molOmega - angularVelocity_[n]) & molMoI)
                                            *(molOmega-angularVelocity_[n]) 
                                            );

                    virialTensor[n] += 0.5*molI->rf();
                }
            }
        }
    }

    //- parallel communication
    if(Pstream::parRun())
    {
        forAll(mass, i)
        {
            forAll(mass[i], j)
            {
                reduce(dof[i][j], sumOp<scalar>());
                reduce(kE[i][j], sumOp<scalar>());
                reduce(angularKeSum[i][j], sumOp<scalar>());
                reduce(kineticTensor[i][j], sumOp<tensor>());
                reduce(virialTensor[i][j], sumOp<tensor>());
            }
        }        
    }


    // collect and compute properties
    

    const scalar& kB = molCloud_.redUnits().kB();

    const scalar& nAvTimeSteps = 1;
    
    forAll(mols, n)
    {
        mols_[n].append(mols[n]);
        mass_[n].append(mass[n]);
        mom_[n].append(mom[n]);
        kE_[n].append(kE[n]);
        angularKeSum_[n].append(angularKeSum[n]);
        dof_[n].append(dof[n]);
        kineticTensor_[n].append(kineticTensor[n]);
        virialTensor_[n].append(virialTensor[n]);
        
        scalar volume = binModel_->binVolume(n);

        if(overideVolume_)
        {
            volume = volumes_[n];
        }
        
        N_[n].append(mols[n]/nAvTimeSteps);
        rhoN_[n].append(mols[n]/(nAvTimeSteps*volume));
        rhoM_[n].append(mass[n]/(nAvTimeSteps*volume));

        if(mass[n] > 0.0)
        {
            UCAM_[n].append(mom[n]/mass[n];
        }
        else
        {
            UCAM_[n].append(vector::zero); 
        }

        if(dof[n] > 0.0)
        {
            T_[n].append((2.0*(kE[n]+angularKeSum[n]))/(kB*dof[n]));
        }
        else
        {
            T_[n].append(0.0);
        }
            

        if(dof[n] > 0.0)
        {
            p_[n].append( tr( (3.0*mols[n]*kineticTensor[n]/dof[n]) + virialTensor[n])
                                    /( 3.0*volume*nAvTimeSteps ) );

            stress_[n].append(((3.0*mols[n]*kineticTensor[n]/dof[n]) + virialTensor[n])
                                    /( volume*nAvTimeSteps ) );
        }
        else
        {
            p_[n].append(0.0);
            stress_[n].append(tensor::zero);
        }
    }
    
}


void polyInstantBinsMethod::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            const reducedUnits& rU = molCloud_.redUnits();            
            scalarField bins = binModel_->binPositions();
            vectorField vectorBins = binModel_->bins();

            label nBins = bins.size();
            label nTimeSteps = N_[0].size();
            
            for (int j = 0; j < nTimeSteps; j++)
            {
                scalarField N(nBins, 0.0);
                scalarField rhoN(nBins, 0.0);
                scalarField rhoM(nBins, 0.0);
                scalarField UCAM_X(nBins, 0.0);
                scalarField UCAM_Y(nBins, 0.0);
                scalarField UCAM_Z(nBins, 0.0);
                scalarField T(nBins, 0.0);
                scalarField p(nBins, 0.0);  
                
                forAll(N_, i)
                {
                    N[i] = N_[i][j];
                    rhoN[i] = rhoN_[i][j];
                    rhoM[i] = rhoM_[i][j];
                    UCAM_X[i] = UCAM_[i][j].x();
                    UCAM_Y[i] = UCAM_[i][j].y();
                    UCAM_Z[i] = UCAM_[i][j].z();
                    T[i] = T_[i][j];
                    p[i] = p_[i][j];
                }
            
                writeTimeData
                (
                    casePath_,
                    "bins_instant_oneDim_"+fieldName_+"_N_SI.xy",
                    N,
                    "sidewaysAppend"
                );

                writeTimeData
                (
                    casePath_,
                    "bins_instant_oneDim_"+fieldName_+"_rhoN_SI.xy",
                    rhoN*rU.refNumberDensity(),
                    "sidewaysAppend"
                );
                
                writeTimeData
                (
                    casePath_,
                    "bins_instant_oneDim_"+fieldName_+"_rhoM_SI.xy",
                    rhoM*rU.refMassDensity(),
                    "sidewaysAppend"
                );
                
                writeTimeData
                (
                    casePath_,
                    "bins_instant_oneDim_"+fieldName_+"_UCAM_X_SI.xy",
                    UCAM_X*rU.refVelocity(),
                    "sidewaysAppend"
                );                    
                writeTimeData
                (
                    casePath_,
                    "bins_instant_oneDim_"+fieldName_+"_UCAM_Y_SI.xy",
                    UCAM_Y*rU.refVelocity(),
                    "sidewaysAppend"
                ); 
                writeTimeData
                (
                    casePath_,
                    "bins_instant_oneDim_"+fieldName_+"_UCAM_Z_SI.xy",
                    UCAM_Z*rU.refVelocity(),
                    "sidewaysAppend"
                );
                
                writeTimeData
                (
                    casePath_,
                    "bins_instant_oneDim_"+fieldName_+"_T_SI.xy",
                    T*rU.refTemp(),
                    "sidewaysAppend"
                );
                
                writeTimeData
                (
                    casePath_,
                    "bins_instant_oneDim_"+fieldName_+"_p_SI.xy",
                    p*rU.refPressure(),
                    "sidewaysAppend"
                );                
            }
        }
    }
}

void polyInstantBinsMethod::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyInstantBinsMethod::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& polyInstantBinsMethod::fields() const
{
    return  fields_;
}

} // End namespace Foam

// ************************************************************************* //
