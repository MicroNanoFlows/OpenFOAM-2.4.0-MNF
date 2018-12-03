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

Measures overall temperature, including vibrational temperature, for a single 
species gas or a gas mixture and writes the results to a volume scalar field 
that can be viewed in Paraview.

Translational, rotatational and vibrational temperature field will also be 
written automatically.

Boundary fields are measured in conjunction with the boundaryMeasurements
class and are also written.

\*---------------------------------------------------------------------------*/

#include "dsmcVolFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcVolFields, 0);

addToRunTimeSelectionTable(dsmcField, dsmcVolFields, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcVolFields::dsmcVolFields
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    n_(),
    t1_(),
    t2_(),
    sampleInterval_(1),
    sampleCounter_(0),
    mfpReferenceTemperature_(273.0),
    fieldName_(propsDict_.lookup("fieldName")),
    dsmcRhoN_
    (
        IOobject
        (
            "dsmcRhoN_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless/dimVolume, 0.0)
    ),
    dsmcRhoNMean_
    (
        IOobject
        (
            "dsmcRhoNMean_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless/dimVolume, 0.0)
    ),
    rhoN_
    (
        IOobject
        (
            "rhoN_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless/dimVolume, 0.0)
    ),
    rhoM_
    (
        IOobject
        (
            "rhoM_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless/dimVolume, 0.0)
    ),
    p_
    (
        IOobject
        (
            "p_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimPressure, 0.0)
    ),
    translationalT_
    (
        IOobject
        (
            "translationalT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    rotationalT_
    (
        IOobject
        (
            "rotationalT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    vibrationalT_
    (
        IOobject
        (
            "vibrationalT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    electronicT_
    (
        IOobject
        (
            "electronicT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    overallT_
    (
        IOobject
        (
            "overallT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    q_
    (
        IOobject
        (
            "surfaceHeatTransfer_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, 0, -3, 0, 0), 0.0)
    ),
    tau_
    (
        IOobject
        (
            "surfaceShearStress_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimPressure, 0.0)
    ),
    MFP_
    (
        IOobject
        (
            "variableHardSphereMeanFreePath_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    MFPdX_
    (
        IOobject
        (
            "meanFreePathCellRatio_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    MCR_
    (
        IOobject
        (
            "meanCollisionRate_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 0, -1, 0, 0), 0.0)
    ),
    MCT_
    (
        IOobject
        (
            "meanCollisionTime_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 0, 1, 0, 0), 0.0)
    ),
    MCTdt_
    (
        IOobject
        (
            "mctTimeStepRatio_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimless, 0.0)
    ),
    MCS_
    (
        IOobject
        (
            "meanCollisionSeparation_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    SOF_
    (
        IOobject
        (
            "separationOfFreePaths_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    Ma_
    (
        IOobject
        (
            "Ma_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    cIDist_
    (
        IOobject
        (
            "classIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    cIIDist_
    (
        IOobject
        (
            "classIIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    cIIIDist_
    (
        IOobject
        (
            "classIIIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    densityError_
    (
        IOobject
        (
            "densityError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    velocityError_
    (
        IOobject
        (
            "velocityError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    temperatureError_
    (
        IOobject
        (
            "temperatureError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    pressureError_
    (
        IOobject
        (
            "pressureError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    UMean_
    (
        IOobject
        (
            "UMean_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0.0", dimLength/dimTime, vector::zero)
    ),
    fD_
    (
        IOobject
        (
            "fD_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0),
            vector::zero
        )
    ),
    heatFluxVector_
    (
        IOobject
        (
            "heatFluxVector_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, 0, -3, 0, 0),
            vector::zero
        )
    ),
    pressureTensor_
    (
        IOobject
        (
            "pressureTensor_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor
        (
            "zero",
            dimPressure,
            tensor::zero
        )
    ),
    shearStressTensor_
    (
        IOobject
        (
            "shearStressTensor_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor
        (
            "zero",
            dimPressure,
            tensor::zero
        )
    ),
    nTimeSteps_(0.0),
    typeIds_(),
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoNInstantaneous_(mesh_.nCells(), 0.0),
    rhoNMeanXnParticle_(mesh_.nCells(), 0.0),
    rhoNMeanInt_(mesh_.nCells(), 0.0),
    molsElec_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    rhoMMeanXnParticle_(mesh_.nCells(), 0.0),
    linearKEMean_(mesh_.nCells(), 0.0),
    linearKEMeanXnParticle_(mesh_.nCells(), 0.0),
    rotationalEMean_(mesh_.nCells(), 0.0),
    rotationalDofMean_(mesh_.nCells(), 0.0),
    muu_(mesh_.nCells(), 0.0),
    muv_(mesh_.nCells(), 0.0),
    muw_(mesh_.nCells(), 0.0),
    mvv_(mesh_.nCells(), 0.0),
    mvw_(mesh_.nCells(), 0.0),
    mww_(mesh_.nCells(), 0.0),
    mcc_(mesh_.nCells(), 0.0),
    mccu_(mesh_.nCells(), 0.0),
    mccv_(mesh_.nCells(), 0.0),
    mccw_(mesh_.nCells(), 0.0),
    eu_(mesh_.nCells(), 0.0),
    ev_(mesh_.nCells(), 0.0),
    ew_(mesh_.nCells(), 0.0),
    e_(mesh_.nCells(), 0.0),
    totalvDof_(mesh_.nCells(), 0.0),
    nClassI_(mesh_.nCells(), 0.0),
    nClassII_(mesh_.nCells(), 0.0),
    nClassIII_(mesh_.nCells(), 0.0),
    collisionSeparation_(mesh_.nCells(), 0.0),
    nColls_(mesh_.nCells(), 0.0),
    momentumMean_(mesh.nCells(), vector::zero),
    momentumMeanXnParticle_(mesh.nCells(), vector::zero),
    boundaryCells_(),
    vibrationalETotal_(),
    electronicETotal_(),
    nParcels_(),
    nParcelsXnParticle_(),
    mccSpecies_(),
    vibT_(),
    nGroundElectronicLevel_(),
    nFirstElectronicLevel_(),
    vDof_(),
    mfp_(),
    mcr_(),
    rhoNBF_(),
    rhoMBF_(),
    linearKEBF_(),
    rotationalEBF_(),
    rotationalDofBF_(),
    qBF_(),
    vibTxvDofBF_(),
    totalvDofBF_(),
    speciesRhoNIntBF_(),
    speciesRhoNElecBF_(),
    momentumBF_(),
    fDBF_(),
    vibrationalEBF_(),
    electronicEBF_(),
    speciesRhoNBF_(),
    mccSpeciesBF_(),
    vibTBF_(),
    vDofBF_(),
    averagingAcrossManyRuns_(false),
    measureClassifications_(false),
    measureMeanFreePath_(false),
    measureErrors_(false),
    densityOnly_(false),
    measureHeatFluxShearStress_(false)
{

    // standard to reading typeIds ------------ 
    const List<word>& molecules (propsDict_.lookup("typeIds"));

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

        if(findIndex(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcVolFields::dsmcVolFields()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << t.system()/"fieldPropertiesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }

    // ---------------------------------------------------

    // Note; outer list is typeIds, inner list is number of cells on the 
    // mesh

    vibT_.setSize(typeIds_.size());

    forAll(vibT_, i)
    {
        vibT_[i].setSize(mesh_.nCells());
    }

    nGroundElectronicLevel_.setSize(typeIds_.size());

    forAll(nGroundElectronicLevel_, i)
    {
        nGroundElectronicLevel_[i].setSize(mesh_.nCells(), 0.0);
    }

    nFirstElectronicLevel_.setSize(typeIds_.size());

    forAll(nFirstElectronicLevel_, i)
    {
        nFirstElectronicLevel_[i].setSize(mesh_.nCells(), 0.0);
    }

    vDof_.setSize(typeIds_.size());

    forAll(vDof_, i)
    {
        vDof_[i].setSize(mesh_.nCells());
    }

    vibrationalETotal_.setSize(typeIds_.size());

    electronicETotal_.setSize(typeIds_.size());

    forAll(electronicETotal_, i)
    {
        electronicETotal_[i].setSize(mesh_.nCells(), 0.0);
    }

    nParcels_.setSize(typeIds_.size());

    forAll(nParcels_, i)
    {
        nParcels_[i].setSize(mesh_.nCells());
    }

    nParcelsXnParticle_.setSize(typeIds_.size());

    forAll(nParcelsXnParticle_, i)
    {
        nParcelsXnParticle_[i].setSize(mesh_.nCells());
    }

    mccSpecies_.setSize(typeIds_.size());

    forAll(nParcels_, i)
    {
        mccSpecies_[i].setSize(mesh_.nCells());
    }

    mfp_.setSize(typeIds_.size());

    forAll(mfp_, i)
    {
        mfp_[i].setSize(mesh_.nCells());
    }

    mcr_.setSize(typeIds_.size());

    forAll(mcr_, i)
    {
        mcr_[i].setSize(mesh_.nCells());
    }

    boundaryCells_.setSize(mesh.boundaryMesh().size());

    forAll(boundaryCells_, p)
    {
        const polyPatch& patch = mesh.boundaryMesh()[p];

        boundaryCells_[p].setSize(patch.size());

        forAll(boundaryCells_[p], c)
        {
            boundaryCells_[p][c] = patch.faceCells()[c];
        }
    }

    // initialisation
    rhoNBF_.setSize(mesh_.boundaryMesh().size());
    rhoMBF_.setSize(mesh_.boundaryMesh().size());
    linearKEBF_.setSize(mesh_.boundaryMesh().size());
    momentumBF_.setSize(mesh_.boundaryMesh().size());
    rotationalEBF_.setSize(mesh_.boundaryMesh().size());
    rotationalDofBF_.setSize(mesh_.boundaryMesh().size());
    qBF_.setSize(mesh_.boundaryMesh().size());
    fDBF_.setSize(mesh_.boundaryMesh().size());
    vibTxvDofBF_.setSize(mesh_.boundaryMesh().size());
    totalvDofBF_.setSize(mesh_.boundaryMesh().size());
    speciesRhoNIntBF_.setSize(mesh_.boundaryMesh().size());
    speciesRhoNElecBF_.setSize(mesh_.boundaryMesh().size());
    n_.setSize(mesh_.boundaryMesh().size());
    t1_.setSize(mesh_.boundaryMesh().size());
    t2_.setSize(mesh_.boundaryMesh().size());

    forAll(rhoNBF_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];

        rhoNBF_[j].setSize(patch.size(), 0.0);
        rhoMBF_[j].setSize(patch.size(), 0.0);
        linearKEBF_[j].setSize(patch.size(), 0.0);
        momentumBF_[j].setSize(patch.size(), vector::zero);
        rotationalEBF_[j].setSize(patch.size(), 0.0);
        rotationalDofBF_[j].setSize(patch.size(), 0.0);
        qBF_[j].setSize(patch.size(), 0.0);
        fDBF_[j].setSize(patch.size(), vector::zero);
        vibTxvDofBF_[j].setSize(patch.size(), 0.0);
        totalvDofBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNIntBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNElecBF_[j].setSize(patch.size(), 0.0);
        n_[j].setSize(patch.size(), vector::zero);
        t1_[j].setSize(patch.size(), vector::zero);
        t2_[j].setSize(patch.size(), vector::zero);
    }

    calculateWallUnitVectors();

    vibrationalEBF_.setSize(typeIds_.size());
    electronicEBF_.setSize(typeIds_.size());
    speciesRhoNBF_.setSize(typeIds_.size());
    mccSpeciesBF_.setSize(typeIds_.size());
    vibTBF_.setSize(typeIds_.size());
    vDofBF_.setSize(typeIds_.size());

    forAll(vibrationalEBF_, i)
    {
        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        electronicEBF_[i].setSize(mesh_.boundaryMesh().size());
        speciesRhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        mccSpeciesBF_[i].setSize(mesh_.boundaryMesh().size());
        vibTBF_[i].setSize(mesh_.boundaryMesh().size());
        vDofBF_[i].setSize(mesh_.boundaryMesh().size());

        forAll(vibrationalEBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            vibrationalEBF_[i][j].setSize(patch.size(), 0.0);
            electronicEBF_[i][j].setSize(patch.size(), 0.0);
            speciesRhoNBF_[i][j].setSize(patch.size(), 0.0);
            mccSpeciesBF_[i][j].setSize(patch.size(), 0.0);
            vibTBF_[i][j].setSize(patch.size(), 0.0);
            vDofBF_[i][j].setSize(patch.size(), 0.0);
        }
    }

    if (propsDict_.found("sampleInterval"))
    {
        sampleInterval_ = 
                        readLabel(propsDict_.lookup("sampleInterval"));
    }
    
    if (propsDict_.found("measureClassifications"))
    {
        measureClassifications_ = 
            Switch(propsDict_.lookup("measureClassifications"));

        if(measureClassifications_)
        {
            Info << "measureClassifications initiated" << endl;
        }
    }

    if (propsDict_.found("measureErrors"))
    {
        measureErrors_ = Switch(propsDict_.lookup("measureErrors"));

        if(measureErrors_)
        {
            Info << "measureErrors initiated" << endl;
        }
    }

    if (propsDict_.found("densityOnly"))
    {
        densityOnly_ = Switch(propsDict_.lookup("densityOnly"));

        if(densityOnly_)
        {
            Info << "densityOnly initiated" << endl;
        }
    }

    if (propsDict_.found("measureHeatFluxShearStress"))
    {
        measureHeatFluxShearStress_ = 
            Switch(propsDict_.lookup("measureHeatFluxShearStress"));

        if(measureHeatFluxShearStress_)
        {
            Info << "measureHeatFluxShearStress initiated" << endl;
        }
    }

    if(propsDict_.found("measureMeanFreePath"))
    {
        measureMeanFreePath_ = Switch(propsDict_.lookup("measureMeanFreePath"));

        if(measureMeanFreePath_)
        {
            Info << "measureMeanFreePath initiated" << endl;

            mfpReferenceTemperature_ = 
            readScalar(propsDict_.lookup("mfpReferenceTemperature"));
        }
    }

    if (propsDict_.found("averagingAcrossManyRuns"))
    {
        averagingAcrossManyRuns_ = 
            Switch(propsDict_.lookup("averagingAcrossManyRuns"));

        // read in stored data from dictionary
        if(averagingAcrossManyRuns_)
        {
            Info << "averagingAcrossManyRuns initiated." << nl << endl;
            readIn();
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcVolFields::~dsmcVolFields()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcVolFields::readIn()
{           
    IOdictionary volFieldsStorage
    (
        IOobject
        (
            "volFieldsMethod_"+fieldName_,
            time_.time().timeName(),
            "uniform",
            time_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    if(volFieldsStorage.size() > 0)
    {
        dictionary
                dict(volFieldsStorage.readStream(volFieldsStorage.filePath()));

        dict.readIfPresent("nTimeSteps", nTimeSteps_);
        dict.readIfPresent("rhoNMean", rhoNMean_);
        dict.readIfPresent("rhoNInstantaneous", rhoNInstantaneous_);
        dict.readIfPresent("rhoNMeanXnParticle", rhoNMeanXnParticle_);
        dict.readIfPresent("rhoNMeanInt", rhoNMeanInt_);
        dict.readIfPresent("molsElec", molsElec_);
        dict.readIfPresent("rhoMMean", rhoMMean_);
        dict.readIfPresent("rhoMMeanXnParticle", rhoMMeanXnParticle_);
        dict.readIfPresent("linearKEMean", linearKEMean_);
        dict.readIfPresent("linearKEMeanXnParticle", linearKEMeanXnParticle_);
        dict.readIfPresent("rotationalEMean", rotationalEMean_);
        dict.readIfPresent("rotationalDofMean", rotationalDofMean_);
        dict.readIfPresent("muu", muu_);
        dict.readIfPresent("muv", muv_);
        dict.readIfPresent("muw", muw_);
        dict.readIfPresent("mvv", mvv_);
        dict.readIfPresent("mvw", mvw_);
        dict.readIfPresent("mww", mww_);
        dict.readIfPresent("mcc", mcc_);
        dict.readIfPresent("mccu", mccu_);
        dict.readIfPresent("mccv", mccv_);
        dict.readIfPresent("mccw", mccw_);
        dict.readIfPresent("eu", eu_);
        dict.readIfPresent("ev", ev_);
        dict.readIfPresent("ew", ew_);
        dict.readIfPresent("e", e_);
        dict.readIfPresent("totalvDof", totalvDof_);
        dict.readIfPresent("nClassI", nClassI_);
        dict.readIfPresent("nClassII", nClassII_);
        dict.readIfPresent("nClassIII", nClassIII_);
        dict.readIfPresent("collisionSeparation", collisionSeparation_);
        dict.readIfPresent("nColls", nColls_);
        dict.readIfPresent("momentumMean", momentumMean_);
        dict.readIfPresent("momentumMeanXnParticle", momentumMeanXnParticle_);
        dict.readIfPresent("vibrationalETotal", vibrationalETotal_);
        dict.readIfPresent("electronicETotal", electronicETotal_);
        dict.readIfPresent("nParcels", nParcels_); 
        dict.readIfPresent("nParcelsXnParticle", nParcelsXnParticle_);
        dict.readIfPresent("mccSpecies", mccSpecies_);
        dict.readIfPresent("vibT", vibT_);
        dict.readIfPresent("nGroundElectronicLevel", nGroundElectronicLevel_);
        dict.readIfPresent("nFirstElectronicLevel", nFirstElectronicLevel_);
        dict.readIfPresent("vDof", vDof_);
        dict.readIfPresent("mfp", mfp_);
        dict.readIfPresent("mcr", mcr_);
        dict.readIfPresent("rhoNBF", rhoNBF_);
        dict.readIfPresent("rhoMBF", rhoMBF_);
        dict.readIfPresent("linearKEBF", linearKEBF_);
        dict.readIfPresent("rotationalEBF", rotationalEBF_);
        dict.readIfPresent("rotationalDofBF", rotationalDofBF_);
        dict.readIfPresent("qBF", qBF_);
        dict.readIfPresent("vibTxvDofBF", vibTxvDofBF_);
        dict.readIfPresent("totalvDofBF", totalvDofBF_);
        dict.readIfPresent("speciesRhoNIntBF", speciesRhoNIntBF_);
        dict.readIfPresent("speciesRhoNElecBF", speciesRhoNElecBF_);
        dict.readIfPresent("momentumBF", momentumBF_);
        dict.readIfPresent("fDBF", fDBF_);
        dict.readIfPresent("vibrationalEBF", vibrationalEBF_);
        dict.readIfPresent("electronicEBF", electronicEBF_);
        dict.readIfPresent("speciesRhoNBF", speciesRhoNBF_);
        dict.readIfPresent("mccSpeciesBF", mccSpeciesBF_);
        dict.readIfPresent("vibTBF", vibTBF_);
        dict.readIfPresent("vDofBF", vDofBF_);  
    }
}

void dsmcVolFields::writeOut()
{
    if (time_.time().outputTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "volFieldsMethod_"+fieldName_,
                time_.time().timeName(),
                "uniform",
                time_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        dict.add("nTimeSteps", nTimeSteps_);
        dict.add("rhoNMean", rhoNMean_);
        dict.add("rhoNInstantaneous", rhoNInstantaneous_);
        dict.add("rhoNMeanXnParticle", rhoNMeanXnParticle_);
        dict.add("rhoNMeanInt", rhoNMeanInt_);
        dict.add("molsElec", molsElec_);
        dict.add("rhoMMean", rhoMMean_);
        dict.add("rhoMMeanXnParticle", rhoMMeanXnParticle_);
        dict.add("linearKEMean", linearKEMean_);
        dict.add("linearKEMeanXnParticle", linearKEMeanXnParticle_);
        dict.add("rotationalEMean", rotationalEMean_);
        dict.add("rotationalDofMean", rotationalDofMean_);
        dict.add("muu", muu_);
        dict.add("muv", muv_);
        dict.add("muw", muw_);
        dict.add("mvv", mvv_);
        dict.add("mvw", mvw_);
        dict.add("mww", mww_);
        dict.add("mcc", mcc_);
        dict.add("mccu", mccu_);
        dict.add("mccv", mccv_);
        dict.add("mccw", mccw_);
        dict.add("eu", eu_);
        dict.add("ev", ev_);
        dict.add("ew", ew_);
        dict.add("e", e_);
        dict.add("totalvDof", totalvDof_);
        dict.add("nClassI", nClassI_);
        dict.add("nClassII", nClassII_);
        dict.add("nClassIII", nClassIII_);
        dict.add("collisionSeparation", collisionSeparation_);
        dict.add("nColls", nColls_);
        dict.add("momentumMean", momentumMean_);
        dict.add("momentumMeanXnParticle", momentumMeanXnParticle_);
        dict.add("vibrationalETotal", vibrationalETotal_);
        dict.add("electronicETotal", electronicETotal_);
        dict.add("nParcels", nParcels_); 
        dict.add("nParcelsXnParticle", nParcelsXnParticle_);
        dict.add("mccSpecies", mccSpecies_);
        dict.add("vibT", vibT_);
        dict.add("nGroundElectronicLevel", nGroundElectronicLevel_);
        dict.add("nFirstElectronicLevel", nFirstElectronicLevel_);
        dict.add("vDof", vDof_);
        dict.add("mfp", mfp_);
        dict.add("mcr", mcr_);
        dict.add("rhoNBF", rhoNBF_);
        dict.add("rhoMBF", rhoMBF_);
        dict.add("linearKEBF", linearKEBF_);
        dict.add("rotationalEBF", rotationalEBF_);
        dict.add("rotationalDofBF", rotationalDofBF_);
        dict.add("qBF", qBF_);
        dict.add("vibTxvDofBF", vibTxvDofBF_);
        dict.add("totalvDofBF", totalvDofBF_);
        dict.add("speciesRhoNIntBF", speciesRhoNIntBF_);
        dict.add("speciesRhoNElecBF", speciesRhoNElecBF_);
        dict.add("momentumBF", momentumBF_);
        dict.add("fDBF", fDBF_);
        dict.add("vibrationalEBF", vibrationalEBF_);
        dict.add("electronicEBF", electronicEBF_);
        dict.add("speciesRhoNBF", speciesRhoNBF_);
        dict.add("mccSpeciesBF", mccSpeciesBF_);
        dict.add("vibTBF", vibTBF_);
        dict.add("vDofBF", vDofBF_);

        IOstream::streamFormat fmt = time_.time().writeFormat();
        IOstream::versionNumber ver = time_.time().writeVersion();
        IOstream::compressionType cmp = time_.time().writeCompression();

        dict.regIOobject::writeObject(fmt, ver, cmp);
    }
}

void dsmcVolFields::calculateWallUnitVectors()
{
    forAll(n_, patchi)
    {
        const polyPatch& pPatch = mesh_.boundaryMesh()[patchi];

        if(isA<wallPolyPatch>(pPatch))
        {
            const vectorField& fC = pPatch.faceCentres();

            forAll(n_[patchi], facei)
            {
                n_[patchi][facei] = pPatch.faceAreas()[facei]
                    /mag(pPatch.faceAreas()[facei]);

                //- Wall tangential unit vector. Use the direction between the
                // face centre and the first vertex in the list
                t1_[patchi][facei] = fC[facei] - 
                            mesh_.points()[mesh_.faces()[pPatch.start() 
                            + facei][0]];
                t1_[patchi][facei] /= mag(t1_[patchi][facei]);

                //- Other tangential unit vector.  Rescaling in case face is not
                //  flat and n and t1 aren't perfectly orthogonal
                t2_[patchi][facei] = n_[patchi][facei]^t1_[patchi][facei];
                t2_[patchi][facei] /= mag(t2_[patchi][facei]);
            }
        }
    }
}

//- initial conditions
void dsmcVolFields::createField()
{
    Info << "Initialising dsmcVolFields field" << endl;

    forAll(vibrationalETotal_, i)
    {
        vibrationalETotal_[i].setSize(cloud_.constProps(typeIds_[i]).
        vibrationalDegreesOfFreedom());

        forAll(vibrationalETotal_[i], j)
        {
            vibrationalETotal_[i][j].setSize(mesh_.nCells(),0.0);
        }
    }
}


void dsmcVolFields::calculateField()
{ 
    sampleCounter_++;

    const scalar& nParticle = cloud_.nParticle();
    rhoNInstantaneous_ = 0.0;

    if(sampleInterval_ <= sampleCounter_)
    {
        nTimeSteps_ += 1.0;

        if(densityOnly_)
        {
            forAllConstIter(dsmcCloud, cloud_, iter)
            {
                const dsmcParcel& p = iter();
                label iD = findIndex(typeIds_, p.typeId());

                if(iD != -1)
                {
                    const label& cell = p.cell();
                    const scalar& mass = cloud_.constProps(p.typeId()).mass();

                    rhoNMean_[cell] += 1.0;
                    rhoNInstantaneous_[cell] += 1.0;

                    if(cloud_.axisymmetric())
                    {
                        const point& cC = cloud_.mesh().cellCentres()[cell];

                        scalar radius = cC.y();
                        scalar RWF = 1.0;
                        RWF = 1.0 + 
                            cloud_.maxRWF()*(radius/cloud_.radialExtent());

                        rhoNMeanXnParticle_[cell] += (RWF*nParticle);
                        rhoMMeanXnParticle_[cell] += (mass*RWF*nParticle);
                    }
                    else
                    {
                        rhoNMeanXnParticle_[cell] += nParticle;
                        rhoMMeanXnParticle_[cell] += (mass*nParticle);
                    }
                }
            }
        }
        else
        {
            forAllConstIter(dsmcCloud, cloud_, iter)
            {
                const dsmcParcel& p = iter();
                const label& iD = findIndex(typeIds_, p.typeId());

                if(iD != -1)
                {
                    const label& cell = p.cell();
                    const scalar& mass = cloud_.constProps(p.typeId()).mass();
                    const scalar& massByMagUsq = mass*sqr(mag(p.U()));
                    const scalarList& electronicEnergies = 
                        cloud_.constProps(typeIds_[iD]).electronicEnergyList();
                    const scalar& rotationalDof = 
                    cloud_.constProps(p.typeId()).rotationalDegreesOfFreedom();
                    const scalar& xVel = p.U().x();
                    const scalar& yVel = p.U().y();
                    const scalar& zVel = p.U().z();
                    const vector& U = p.U();

                    scalarList EVib
                    (
                        cloud_.constProps(typeIds_[iD])
                        .vibrationalDegreesOfFreedom()
                    );

                    if(EVib.size() > 0)
                    {
                        forAll(EVib, i)
                        {

                            EVib[i] = p.vibLevel()[i]
                                     *physicoChemical::k.value()                
                                     *cloud_.constProps(p.typeId()).thetaV()[i];

                            vibrationalETotal_[iD][i][cell] +=
                                        p.vibLevel()[i]
                                        *physicoChemical::k.value()
                                     *cloud_.constProps(p.typeId()).thetaV()[i];
                        }
                    }

                    rhoNMean_[cell] += 1.0;
                    rhoNInstantaneous_[cell] += 1.0;
                    rhoMMean_[cell] += mass;
                    linearKEMean_[cell] += mass*(U & U);
                    momentumMean_[cell] += mass*U;
                    rotationalEMean_[cell] += p.ERot();
                    rotationalDofMean_[cell] += rotationalDof;
                    electronicETotal_[iD][cell] += 
                                        electronicEnergies[p.ELevel()];
                    nParcels_[iD][cell] += 1.0;
                    mccSpecies_[iD][cell] += massByMagUsq;

                    if(cloud_.axisymmetric())
                    {
                        const point& cC = cloud_.mesh().cellCentres()[cell];

                        scalar radius = cC.y();
                        scalar RWF = 1.0;
                        RWF = 1.0 + 
                            cloud_.maxRWF()*(radius/cloud_.radialExtent());

                        nParcelsXnParticle_[iD][cell] += (RWF*nParticle);
                        rhoNMeanXnParticle_[cell] += (RWF*nParticle);
                        rhoMMeanXnParticle_[cell] += (mass*RWF*nParticle);
                        momentumMeanXnParticle_[cell] += 
                                    (mass*(U)*RWF*nParticle);
                        linearKEMeanXnParticle_[cell] += (mass*(U & 
                                                        U)*RWF*nParticle);
                    }
                    else
                    {
                        nParcelsXnParticle_[iD][cell] += nParticle;
                        rhoNMeanXnParticle_[cell] += nParticle;
                        rhoMMeanXnParticle_[cell] += (mass*nParticle);
                        momentumMeanXnParticle_[cell] += 
                                                (mass*(U)*nParticle);
                        linearKEMeanXnParticle_[cell] += (mass*(U & 
                                                        U)*nParticle);
                    }

                    muu_[cell] += mass*sqr(xVel);
                    muv_[cell] += mass*( xVel * yVel );
                    muw_[cell] += mass*( xVel * zVel );
                    mvv_[cell] += mass*sqr(yVel);
                    mvw_[cell] += mass*( yVel * zVel );
                    mww_[cell] += mass*sqr(zVel);
                    mcc_[cell] += massByMagUsq;
                    mccu_[cell] += massByMagUsq*(xVel);
                    mccv_[cell] += massByMagUsq*(yVel);
                    mccw_[cell] += massByMagUsq*(xVel);

                    scalar vibEn = 0.0;

                    if(EVib.size() > 0)
                    {
                        forAll(EVib, v)
                        {
                            vibEn += EVib[v];
                        }
                    }
                    
                    eu_[cell] += ( p.ERot() + vibEn )*(xVel);
                    ev_[cell] += ( p.ERot() + vibEn )*(yVel);
                    ew_[cell] += ( p.ERot() + vibEn )*(zVel);
                    e_[cell] += ( p.ERot() + vibEn );
                    
                    if(rotationalDof > VSMALL)
                    {
                        rhoNMeanInt_[cell] += 1.0;
                    }
                    
                    label nElecLevels =                 
                       cloud_.constProps(p.typeId()).numberOfElectronicLevels();
                    
                    if(nElecLevels > 1)
                    {
                        molsElec_[cell] += 1.0;
                        
                        if(p.ELevel() == 0)
                        {
                            nGroundElectronicLevel_[iD][cell]++;
                        }
                        if(p.ELevel() == 1)
                        {
                            nFirstElectronicLevel_[iD][cell]++;
                        }
                    }
                    
                    if(measureClassifications_)
                    {
                        label classification = p.classification();
                        
                        if(classification == 0)
                        {
                            nClassI_[cell] += 1.0;
                        }
                        
                        if(classification == 1)
                        {
                            nClassII_[cell] += 1.0;
                        }
                        
                        if(classification == 2)
                        {
                            nClassIII_[cell] += 1.0;
                        }
                    }
                }
            }
            
            // obtain collision quality measurements
            
            forAll(cloud_.cellPropMeasurements().collisionSeparation(), cell)
            {
                collisionSeparation_[cell] += 
                    cloud_.cellPropMeasurements().collisionSeparation()[cell];
                nColls_[cell] += cloud_.cellPropMeasurements().nColls()[cell];
            }
            
            // obtain boundary measurements
            
            forAll(cloud_.boundaryFluxMeasurements().rhoNBF(), i)
            {
                const label& iD = findIndex(typeIds_, i);
            
                forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i], j)
                {                
                    forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i][j], k)
                    {
                        if(iD != -1)
                        { 
                            rhoNBF_[j][k] += 
                            cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k];
                            rhoMBF_[j][k] += 
                            cloud_.boundaryFluxMeasurements().rhoMBF()[i][j][k];
                            linearKEBF_[j][k] += 
                        cloud_.boundaryFluxMeasurements().linearKEBF()[i][j][k];
                            momentumBF_[j][k] += 
                        cloud_.boundaryFluxMeasurements().momentumBF()[i][j][k];
                            rotationalEBF_[j][k] += 
                    cloud_.boundaryFluxMeasurements().rotationalEBF()[i][j][k];
                            rotationalDofBF_[j][k] += 
                cloud_.boundaryFluxMeasurements().rotationalDofBF()[i][j][k];
                            qBF_[j][k] += 
                               cloud_.boundaryFluxMeasurements().qBF()[i][j][k];
                            fDBF_[j][k] += 
                              cloud_.boundaryFluxMeasurements().fDBF()[i][j][k];
                            speciesRhoNBF_[iD][j][k] += 
                            cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k];
                            vibrationalEBF_[iD][j][k] += 
                    cloud_.boundaryFluxMeasurements().vibrationalEBF()[i][j][k];
                            electronicEBF_[iD][j][k] += 
                    cloud_.boundaryFluxMeasurements().electronicEBF()[i][j][k];
                            mccSpeciesBF_[iD][j][k] += 
                    cloud_.boundaryFluxMeasurements().mccSpeciesBF()[i][j][k];
                            speciesRhoNIntBF_[j][k] += 
                        cloud_.boundaryFluxMeasurements().rhoNIntBF()[i][j][k];
                            speciesRhoNElecBF_[j][k] += 
                        cloud_.boundaryFluxMeasurements().rhoNElecBF()[i][j][k];
                        }
                    }
                }
            }
        }
        sampleCounter_ = 0;
    }
    
    if(time_.time().outputTime())
    {
        const scalar& nAvTimeSteps = nTimeSteps_;
                
        if(densityOnly_)
        {
            forAll(rhoNMean_, cell)
            {
                if(rhoNMean_[cell] > VSMALL)
                {
                    const scalar& cellVolume = mesh_.cellVolumes()[cell];
                    
                    dsmcRhoNMean_[cell] = rhoNMean_[cell]/(nAvTimeSteps);
                    
                    rhoN_[cell] = 
                        (rhoNMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);
                    
                    rhoM_[cell] = 
                        (rhoMMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);

                }
                else
                {
                     // not zero so that weighted decomposition still works
                    dsmcRhoNMean_[cell] = 0.001;
                    rhoN_[cell] = 0.0;
                    rhoM_[cell] = 0.0;
                }
                
                if(rhoNInstantaneous_[cell] > VSMALL)
                {
                    dsmcRhoN_[cell] = rhoNInstantaneous_[cell];
                }
                else
                {
                    dsmcRhoN_[cell] = 0.001;
                }
            }
        }
        else
        {                  
            scalarField vibT(mesh_.nCells(), scalar(0.0));
            scalarField vibTForOverallT(mesh_.nCells(), scalar(0.0));
            scalarField Cp(mesh_.nCells(), scalar(0.0));
            scalarField Cv(mesh_.nCells(), scalar(0.0));
            scalarField molecularMass(mesh_.nCells(), scalar(0.0));
            scalarField Cv_p(mesh_.nCells(), scalar(0.0));
            scalarField totalvDof(mesh_.nCells(), scalar(0.0));
            scalarField totalvDofOverall(mesh_.nCells(), scalar(0.0));
            
            forAll(rhoNMean_, cell)
            {                
                if(rhoNMean_[cell] > VSMALL)
                {                                     
                    const scalar& cellVolume = mesh_.cellVolumes()[cell];
                    
                    dsmcRhoNMean_[cell] = rhoNMean_[cell]/(nAvTimeSteps);
                    
                    rhoN_[cell] = 
                        (rhoNMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);
                    
                    rhoM_[cell] = 
                        (rhoMMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);
                    
                    scalar rhoMMean = 
                        rhoMMeanXnParticle_[cell]/(cellVolume*nAvTimeSteps);
                    UMean_[cell] = momentumMeanXnParticle_[cell] / 
                                    (rhoMMean*cellVolume*nAvTimeSteps);
                    scalar linearKEMean = 0.5*linearKEMeanXnParticle_[cell] 
                                            / (cellVolume*nAvTimeSteps);
                    scalar rhoNMean = 
                        rhoNMeanXnParticle_[cell]/(cellVolume*nAvTimeSteps);

                    translationalT_[cell] = 
                        2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                                    *(linearKEMean - 0.5*rhoMMean*
                                        (UMean_[cell] & UMean_[cell])
                                    );
                                    
                    p_[cell] = 
                        rhoN_[cell]*physicoChemical::k.value()
                        *translationalT_[cell];
                }
                else
                {
                    // not zero so that weighted decomposition still works
                    dsmcRhoNMean_[cell] = 0.001;
                    rhoN_[cell] = 0.0;
                    rhoM_[cell] = 0.0;
                    UMean_[cell] = vector::zero;
                    translationalT_[cell] = 0.0;           
                    p_[cell] = 0.0;
                }
                
                if(rhoNInstantaneous_[cell] > VSMALL)
                {
                    dsmcRhoN_[cell] = rhoNInstantaneous_[cell];
                }
                else
                {
                    dsmcRhoN_[cell] = 0.001;
                }

                if(rotationalDofMean_[cell] > VSMALL && nAvTimeSteps > VSMALL)
                {
                    scalar rotationalEMean = rotationalEMean_[cell] / 
                                                nAvTimeSteps;
                    scalar rotationalDofMean = rotationalDofMean_[cell] / 
                                                nAvTimeSteps;

                    rotationalT_[cell] =  (2.0/physicoChemical::k.value())
                                           *(rotationalEMean/rotationalDofMean);
                }
                else
                {
                    rotationalT_[cell] = 0.0;
                }
               
                scalarList degreesOfFreedomSpecies(typeIds_.size(),0.0);
                scalarList vibTID(vibrationalETotal_.size(), 0.0);
                
                List<scalarList> dofMode;
                List<scalarList> vibTMode;
                
                dofMode.setSize(typeIds_.size());
                vibTMode.setSize(typeIds_.size());
                
                forAll(dofMode, iD)
                { 
                    dofMode[iD].setSize(cloud_.constProps(typeIds_[iD]).
                                        vibrationalDegreesOfFreedom(), 0.0);
                    
                    vibTMode[iD].setSize(cloud_.constProps(typeIds_[iD]).
                                        vibrationalDegreesOfFreedom(), 0.0);
                }
               
                forAll(vibrationalETotal_, iD)
                {
                    forAll(vibrationalETotal_[iD], v)
                    {
                        if(vibrationalETotal_[iD][v][cell] > VSMALL
                                    && nParcels_[iD][cell] > VSMALL
                                    && dofMode.size() > VSMALL)
                        {        
                            scalar thetaV = 
                                cloud_.constProps(typeIds_[iD]).thetaV()[v];
                            
                            scalar vibrationalEMean = 
                                    vibrationalETotal_[iD][v][cell]
                                    /nParcels_[iD][cell];
                            
                            scalar iMean = 
                                    vibrationalEMean
                                    /(physicoChemical::k.value()*thetaV);
                            
                            vibTMode[iD][v] = thetaV / log(1.0 + (1.0/iMean));

                            dofMode[iD][v] = 
                                (2.0*thetaV/vibTMode[iD][v]) 
                                / (exp(thetaV/vibTMode[iD][v]) - 1.0);
                        }
                    }
                    
                    forAll(dofMode[iD], v)
                    {
                        degreesOfFreedomSpecies[iD] += 
                                    dofMode[iD][v];
                    }
                    
                    forAll(dofMode[iD], v)
                    {
                        if(degreesOfFreedomSpecies[iD] > VSMALL)
                        {
                            vibTID[iD] += 
                                vibTMode[iD][v]
                                *(dofMode[iD][v]
                                /degreesOfFreedomSpecies[iD]);
                        }
                    }

                    
                    totalvDof[cell] += degreesOfFreedomSpecies[iD];
                    
                    if(rhoNMeanInt_[cell] > VSMALL 
                        && rhoNMean_[cell] > VSMALL 
                        && nParcels_[iD][cell] > VSMALL)
                    {
                        scalar fraction = 
                                nParcels_[iD][cell]
                                /rhoNMeanInt_[cell];
                        
                        scalar fractionOverall = 
                                nParcels_[iD][cell]
                                /rhoNMean_[cell];
                        
                        totalvDofOverall[cell] += 
                                totalvDof[cell]
                                *(fractionOverall/fraction);
                        
                        vibT[cell] += vibTID[iD]*fraction;
                    }
                }

                vibrationalT_[cell] = vibT[cell];
                                
                // electronic temperature
                scalar totalEDof = 0.0;
                scalar elecT = 0.0;
                    
                forAll(nParcels_, iD)
                {
                    const scalarList& electronicEnergies = 
                        cloud_.constProps(typeIds_[iD]).electronicEnergyList();
                    const labelList& degeneracies = 
                        cloud_.constProps(typeIds_[iD]).degeneracyList();

                    if(nGroundElectronicLevel_[iD][cell] > VSMALL && 
                        nFirstElectronicLevel_[iD][cell] > VSMALL && 
                        nFirstElectronicLevel_[iD][cell]*degeneracies[0] != 
                        nGroundElectronicLevel_[iD][cell]*degeneracies[1])
                    {
                        
                        scalar elecTID = 
                            (electronicEnergies[1]-electronicEnergies[0])/
                            (
                                physicoChemical::k.value()*
                                log((nGroundElectronicLevel_[iD][cell]*         
                                 degeneracies[1])/
                                (nFirstElectronicLevel_[iD][cell]*
                                degeneracies[0]))   
                            );

                    
                        scalar fraction = nParcels_[iD][cell]/molsElec_[cell];
                            
                        if(elecTID > VSMALL)
                        {
                            elecT += fraction*elecTID;
                        }
                        
                        
                        scalar eDof = (2.0*(electronicETotal_[iD][cell]/
                             nParcels_[iD][cell]))/
                             (physicoChemical::k.value()*elecTID);

                        totalEDof += fraction*eDof;
                    }
                }

                electronicT_[cell] = elecT;

                scalar nRotDof = 0.0;
                
                if(rhoNMean_[cell] > VSMALL)
                {
                    nRotDof = rotationalDofMean_[cell] / rhoNMean_[cell];
                }
                
                overallT_[cell] = ( 
                                    (3.0*translationalT_[cell]) 
                                    + (nRotDof*rotationalT_[cell]) 
                                    + (totalvDof_[cell]*vibrationalT_[cell])
                                    + (totalEDof*electronicT_[cell])
                                   ) /
                                    (3.0 + nRotDof + totalvDof_[cell] + 
                                        totalEDof);                 
                                    
                if(measureHeatFluxShearStress_)
                {                    
                    if(rhoNMean_[cell] > VSMALL)
                    {
                        pressureTensor_[cell].xx() = rhoN_[cell]*
                        ( 
                            muu_[cell]/(rhoNMean_[cell]) - 
                            (
                                (rhoMMean_[cell]/(rhoNMean_[cell]))
                                *UMean_[cell].x()*UMean_[cell].x()
                            )
                        );
                        pressureTensor_[cell].xy() = rhoN_[cell]*
                        ( 
                            muv_[cell]/(rhoNMean_[cell]) - 
                            ((rhoMMean_[cell]/(rhoNMean_[cell])))
                            *UMean_[cell].x()*UMean_[cell].y() 
                            
                        );
                        pressureTensor_[cell].xz() = rhoN_[cell]*
                        ( 
                            muw_[cell]/(rhoNMean_[cell]) - 
                            ((rhoMMean_[cell]/(rhoNMean_[cell]))
                            *UMean_[cell].x()*UMean_[cell].z())  
                        );
                        pressureTensor_[cell].yx() = pressureTensor_[cell].xy();
                        pressureTensor_[cell].yy() = rhoN_[cell]*
                        ( 
                            mvv_[cell]/(rhoNMean_[cell]) - 
                            ((rhoMMean_[cell]/(rhoNMean_[cell])))
                            *UMean_[cell].y()*UMean_[cell].y()
                        );
                        pressureTensor_[cell].yz() = rhoN_[cell]*
                        ( 
                            mvw_[cell]/(rhoNMean_[cell]) - 
                            ((rhoMMean_[cell]/(rhoNMean_[cell]))
                            *UMean_[cell].y()*UMean_[cell].z())
                        );
                        pressureTensor_[cell].zx() = pressureTensor_[cell].xz();
                        pressureTensor_[cell].zy() = pressureTensor_[cell].yz();
                        pressureTensor_[cell].zz() = rhoN_[cell]*
                        (
                            mww_[cell]/(rhoNMean_[cell]) - 
                            ((rhoMMean_[cell]/(rhoNMean_[cell]))
                            *UMean_[cell].z()*UMean_[cell].z()) 
                        );
                                                
                        scalar scalarPressure = (1.0/3.0)*
                                                (pressureTensor_[cell].xx() + 
                                                pressureTensor_[cell].yy() + 
                                                pressureTensor_[cell].zz());
                        
                        shearStressTensor_[cell] = -pressureTensor_[cell];
                        shearStressTensor_[cell].xx() += scalarPressure;
                        shearStressTensor_[cell].yy() += scalarPressure;
                        shearStressTensor_[cell].zz() += scalarPressure;
                        
                       //terms involving pressure tensor should not be 
                       //multiplied by the number density (see Bird corrigendum)
                        
                        heatFluxVector_[cell].x() = rhoN_[cell]*
                        (
                            0.5*(mccu_[cell]/(rhoNMean_[cell])) - 
                            0.5*(mcc_[cell]/(rhoNMean_[cell]))*
                            UMean_[cell].x() + eu_[cell]/(rhoNMean_[cell]) - 
                            (e_[cell]/(rhoNMean_[cell]))*UMean_[cell].x()
                        ) - 
                            pressureTensor_[cell].xx()*UMean_[cell].x() - 
                            pressureTensor_[cell].xy()*UMean_[cell].y() - 
                            pressureTensor_[cell].xz()*UMean_[cell].z();

                        heatFluxVector_[cell].y() = rhoN_[cell]*
                        (
                            0.5*(mccv_[cell]/(rhoNMean_[cell])) - 
                            0.5*(mcc_[cell]/(rhoNMean_[cell]))*
                            UMean_[cell].y() + ev_[cell]/(rhoNMean_[cell])-
                            (e_[cell]/(rhoNMean_[cell]))*UMean_[cell].y()
                        ) - 
                            pressureTensor_[cell].yx()*UMean_[cell].x() -
                            pressureTensor_[cell].yy()*UMean_[cell].y() -
                            pressureTensor_[cell].yz()*UMean_[cell].z();

                        heatFluxVector_[cell].z() = rhoN_[cell]*
                        (
                            0.5*(mccw_[cell]/(rhoNMean_[cell])) -
                            0.5*(mcc_[cell]/(rhoNMean_[cell]))*
                            UMean_[cell].z() + ew_[cell]/(rhoNMean_[cell]) -
                            (e_[cell]/(rhoNMean_[cell]))*UMean_[cell].z()
                        ) -
                            pressureTensor_[cell].zx()*UMean_[cell].x() -
                            pressureTensor_[cell].zy()*UMean_[cell].y() -
                            pressureTensor_[cell].zz()*UMean_[cell].z();
                    }
                    else
                    {
                        pressureTensor_[cell] = tensor::zero;
                        shearStressTensor_[cell] = tensor::zero;
                        heatFluxVector_[cell] = vector::zero;
                    }
                }

                totalvDof_[cell] = 0.0;

                forAll(nParcels_, iD)
                {
                    label typeId = typeIds_[iD];

                    if(rhoNMean_[cell] > VSMALL)
                    {
                        molecularMass[cell] += 
                                        cloud_.constProps(typeId).mass()*
                                        (nParcels_[iD][cell]/rhoNMean_[cell]);
                        Cp[cell] += (5.0 + 
                         cloud_.constProps(typeId).rotationalDegreesOfFreedom())
                            *(nParcels_[iD][cell]/rhoNMean_[cell]);
                        Cv[cell] += (3.0 + 
                        cloud_.constProps(typeId).rotationalDegreesOfFreedom())
                            *(nParcels_[iD][cell]/rhoNMean_[cell]);
                    }
                }

                Cv_p[cell] = Cv[cell]/6.02214129e23;

                scalar gasConstant = 0.0;
                scalar gamma = 0.0;
                scalar speedOfSound = 0.0;

                if(molecularMass[cell] > VSMALL)
                {
                    gasConstant = // R = k/m
                        physicoChemical::k.value()/molecularMass[cell]; 
                }

                if(Cv[cell] > VSMALL)
                {
                    gamma = Cp[cell]/Cv[cell]; // gamma = cP/cV
                }

                if(gamma > VSMALL && gasConstant > VSMALL 
                    && translationalT_[cell] > VSMALL)
                {
                    speedOfSound = 
                        sqrt(gamma*gasConstant*translationalT_[cell]);
                }

                if(speedOfSound > VSMALL)
                {
                    Ma_[cell] = mag(UMean_[cell])/speedOfSound;
                }
                else
                {
                    Ma_[cell] = 0.0;
                }

                if(measureMeanFreePath_)
                {
                    forAll(mfp_, iD)
                    {
                        label qspec = 0;

                        for (qspec=0; qspec<typeIds_.size(); qspec++)
                        {
                            scalar dPQ = 
                                0.5*(cloud_.constProps(typeIds_[iD]).d() + 
                                cloud_.constProps(typeIds_[qspec]).d());

                            scalar omegaPQ = 
                                0.5*(cloud_.constProps(typeIds_[iD]).omega() + 
                                cloud_.constProps(typeIds_[qspec]).omega());

                            scalar massRatio = 
                                cloud_.constProps(typeIds_[iD]).mass()/
                                cloud_.constProps(typeIds_[qspec]).mass();

                            if(nParcels_[qspec][cell] > VSMALL && 
                                        translationalT_[cell] > VSMALL)
                            {
                                scalar nDensQ = 
                                    (nParcelsXnParticle_[qspec][cell])/
                                    (mesh_.cellVolumes()[cell]*nTimeSteps_);

                                scalar reducedMass = 
                                (
                                    cloud_.constProps(typeIds_[iD]).mass()*
                                    cloud_.constProps(typeIds_[qspec]).mass()  
                                )/ 
                                (
                                    cloud_.constProps(typeIds_[iD]).mass()+
                                    cloud_.constProps(typeIds_[qspec]).mass() 
                                );

                                mfp_[iD][cell] += 
                                    (
                                        pi*dPQ*dPQ*nDensQ*
                                        pow(mfpReferenceTemperature_/
                                        translationalT_[cell], omegaPQ-0.5)
                                        *sqrt(1.0+massRatio)
                                    ); //Bird, eq (4.76)

                                mcr_[iD][cell] += 
                                (
                                    2.0*sqrt(pi)*dPQ*dPQ*nDensQ*
                                    pow(translationalT_[cell]/
                                    mfpReferenceTemperature_, 1.0-omegaPQ)  
                                    *sqrt(2.0*physicoChemical::k.value()*
                                    mfpReferenceTemperature_/reducedMass)
                                ); // //Bird, eq (4.74)
                            }
                        }

                        if(mfp_[iD][cell] > VSMALL)
                        {
                            mfp_[iD][cell] = 1.0/mfp_[iD][cell];
                        }
                    }

                    MFP_[cell] = 0.0;
                    MFPdX_[cell] = 0.0;
                    MCR_[cell] = 0.0;
                    MCT_[cell] = 0.0;
                    MCTdt_[cell] = 0.0;
                    MCS_[cell] = 0.0;

                    if(nColls_[cell] > VSMALL)
                    {
                        MCS_[cell] = collisionSeparation_[cell]/nColls_[cell];
                    }
                    else
                    {
                       MCS_[cell] = GREAT; 
                    }

                    forAll(mfp_, iD)
                    {
                        if(rhoN_[cell] > VSMALL)
                        {
                            scalar nDensP = (nParcelsXnParticle_[iD][cell])/
                                        (mesh_.cellVolumes()[cell]*nTimeSteps_);

                            //Bird, eq (4.77)
                            MFP_[cell] += mfp_[iD][cell]*nDensP/rhoN_[cell];
 
                            //Bird, eq (1.38)
                            MCR_[cell] += mcr_[iD][cell]*nDensP/rhoN_[cell];
                        }
                    }

                    if(MFP_[cell] < VSMALL)
                    {
                        MFP_[cell] = GREAT;
                    }

                    const scalar& deltaT = mesh_.time().deltaTValue();

                    if(MCR_[cell] > VSMALL)
                    {
                        MCT_[cell] = 1.0/MCR_[cell];
                        MCTdt_[cell] = MCT_[cell]/deltaT;
                    }
                    else
                    {
                        MCT_[cell] = GREAT;
                        MCTdt_[cell] = GREAT;
                    }

                    forAll(mfp_, iD)
                    {
                        mfp_[iD][cell] = 0.0;
                        mcr_[iD][cell] = 0.0;
                    }

                    if(MFP_[cell] != GREAT)
                    {
                        scalar largestCellDimension = 0.0;

                        const labelList& pLabels(mesh_.cells()[cell].
                                                    labels(mesh_.faces()));
                        pointField pLocal(pLabels.size(), vector::zero);

                        forAll (pLabels, pointi)
                        {
                            pLocal[pointi] = mesh_.points()[pLabels[pointi]];
                        }

                        scalarField dimension;

                        dimension.setSize(2, 0.0);

                        dimension[0] =  Foam::max(pLocal & vector(1,0,0)) - 
                                            Foam::min(pLocal & vector(1,0,0));
                        dimension[1] = Foam::max(pLocal & vector(0,1,0)) - 
                                            Foam::min(pLocal & vector(0,1,0));

                        largestCellDimension = dimension[0];

                        label dim = 0;

                        for (dim=0; dim<dimension.size(); dim++)
                        {
                            if(dimension[dim] > largestCellDimension)
                            {
                                largestCellDimension = dimension[dim];
                            }
                        }

                        MFPdX_[cell] = MFP_[cell]/largestCellDimension;

                        if(MFP_[cell] > VSMALL && MCS_[cell] > VSMALL)
                        {
                            SOF_[cell] = MCS_[cell]/MFP_[cell];
                        }
                    }
                    else
                    {
                        MFPdX_[cell] = GREAT;
                        SOF_[cell] = GREAT;
                    }
                }

                if(measureClassifications_)
                {
                    if(rhoNMean_[cell] > VSMALL)
                    {
                        cIDist_[cell] = nClassI_[cell]/rhoNMean_[cell];
                        cIIDist_[cell] = nClassII_[cell]/rhoNMean_[cell];
                        cIIIDist_[cell] = nClassIII_[cell]/rhoNMean_[cell];
                    }
                }

                if(measureErrors_)
                {
                    if(dsmcRhoNMean_[cell] > VSMALL &&
                        Ma_[cell] > VSMALL && gamma > VSMALL
                        && Cv_p[cell] > VSMALL)
                    {
                        densityError_[cell] =1.0/
                                        sqrt(dsmcRhoNMean_[cell]*nTimeSteps_);
                        velocityError_[cell] =
                            (1.0/sqrt(dsmcRhoNMean_[cell]*nTimeSteps_))
                            *(1.0/(Ma_[cell]*sqrt(gamma)));
                        temperatureError_[cell] = 
                                (1.0/sqrt(dsmcRhoNMean_[cell]*nTimeSteps_))*
                                    sqrt(physicoChemical::k.value()/Cv_p[cell]);
                        pressureError_[cell] = sqrt(gamma)/
                                        sqrt(dsmcRhoNMean_[cell]*nTimeSteps_);
                    }
                }
            }

            List<scalarField> vibTBF(mesh_.boundaryMesh().size());
            List<scalarField> molMassBoundary(mesh_.boundaryMesh().size());
            List<scalarField> CpBoundary(mesh_.boundaryMesh().size());
            List<scalarField> CvBoundary(mesh_.boundaryMesh().size());
            List<scalarField> Cv_pBoundary(mesh_.boundaryMesh().size());

            // computing boundary measurements
            forAll(rhoNBF_, j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];

                vibTBF[j].setSize(patch.size(), 0.0);
                molMassBoundary[j].setSize(patch.size(), 0.0);
                CpBoundary[j].setSize(patch.size(), 0.0);
                CvBoundary[j].setSize(patch.size(), 0.0); 
                Cv_pBoundary[j].setSize(patch.size(), 0.0);

                if(isA<wallPolyPatch>(patch))
                {
                    forAll(rhoN_.boundaryField()[j], k)
                    {
                        rhoN_.boundaryField()[j][k] =
                                rhoNBF_[j][k]*nParticle/nAvTimeSteps;
                        rhoM_.boundaryField()[j][k] =
                                rhoMBF_[j][k]*nParticle/nAvTimeSteps;

                        if(rhoM_.boundaryField()[j][k] > VSMALL)
                        {
                            UMean_.boundaryField()[j][k] =
                                    momentumBF_[j][k]*nParticle/
                                    (rhoM_.boundaryField()[j][k]*nAvTimeSteps);
                        }
                        else
                        {
                            UMean_.boundaryField()[j][k] = vector::zero;
                        }

                        scalar rhoMMean = rhoMBF_[j][k]*nParticle/nAvTimeSteps;
                        scalar linearKEMean = linearKEBF_[j][k]*
                                                nParticle/nAvTimeSteps;
                        scalar rhoNMean = rhoNBF_[j][k]*nParticle/nAvTimeSteps;

                        if(rhoNMean > VSMALL)
                        {
                            translationalT_.boundaryField()[j][k] = 
                                2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                                            *(linearKEMean - 0.5*rhoMMean*
                                            (
                                                UMean_.boundaryField()[j][k] & 
                                                UMean_.boundaryField()[j][k])
                                            );
                        }
                        else
                        {
                            translationalT_.boundaryField()[j][k] = 0.0;
                        }

                        if(rotationalDofBF_[j][k] > VSMALL)
                        {
                            rotationalT_.boundaryField()[j][k] = 
                                        (2.0/physicoChemical::k.value())
                                *(rotationalEBF_[j][k]/rotationalDofBF_[j][k]);
                        }
                        else
                        {
                            rotationalT_.boundaryField()[j][k] = 0.0;
                        }

                        // electronic temperature
                        scalar totalEDof = 0.0;
                        scalar elecT = 0.0;

                        electronicT_.boundaryField()[j][k] = elecT;

                        scalar nRotDof = 0.0;

                        if(rhoNBF_[j][k] > VSMALL)
                        {
                            nRotDof = rotationalDofBF_[j][k] / rhoNBF_[j][k];
                        }

                        overallT_.boundaryField()[j][k] = 
                        (
                            (3.0*translationalT_.boundaryField()[j][k]) +
                            (nRotDof*rotationalT_.boundaryField()[j][k]) +
                                (totalvDofBF_[j][k]*
                                vibrationalT_.boundaryField()[j][k]) +
                            (totalEDof*elecT)
                        ) /
                        (3.0 + nRotDof + totalvDofBF_[j][k]+ totalEDof);

                        totalvDofBF_[j][k] = 0.0;

                        Cv_pBoundary[j][k] = CvBoundary[j][k]/6.02214129e23;

                        scalar gasConstant = 0.0;
                        scalar gamma = 0.0;
                        scalar speedOfSound = 0.0;

                        if(molMassBoundary[j][k] > VSMALL)
                        {
                            gasConstant =   // R = k/m
                               physicoChemical::k.value()/molMassBoundary[j][k];
                        }

                        if(CvBoundary[j][k] > VSMALL)
                        {
                            gamma = CpBoundary[j][k]/CvBoundary[j][k];
                        }

                        if(gamma > VSMALL && gasConstant > VSMALL && 
                            translationalT_.boundaryField()[j][k] > VSMALL)
                        {
                            speedOfSound = sqrt
                                    (
                                            gamma*gasConstant*                   
                                        translationalT_.boundaryField()[j][k]
                                    );
                        }

                        if(speedOfSound > VSMALL)
                        {
                            Ma_.boundaryField()[j][k] = 
                                        mag(UMean_.boundaryField()[j][k])/
                                        speedOfSound;
                        }
                        else
                        {
                            Ma_.boundaryField()[j][k] = 0.0;
                        }

                        q_.boundaryField()[j][k] = qBF_[j][k]/nAvTimeSteps;
                        fD_.boundaryField()[j][k] = fDBF_[j][k]/nAvTimeSteps;
                    }

                    p_.boundaryField()[j] = fD_.boundaryField()[j] & n_[j];

                    tau_.boundaryField()[j] = sqrt(
                        sqr(fD_.boundaryField()[j] & t1_[j])
                        + sqr(fD_.boundaryField()[j] & t2_[j]));
                }
            }

            forAll(boundaryCells_, j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];

                labelList bCs = boundaryCells_[j];

                forAll(bCs, k)
                {
                    if(isA<polyPatch>(patch))
                    {
                        if(!isA<emptyPolyPatch>(patch))
                        {
                            if(!isA<cyclicPolyPatch>(patch))
                            {
                                dsmcRhoNMean_.boundaryField()[j][k] = 
                                                    dsmcRhoNMean_[bCs[k]];
                                rhoN_.boundaryField()[j][k] = 
                                                    rhoN_[bCs[k]];
                                rhoM_.boundaryField()[j][k] = 
                                                    rhoM_[bCs[k]];

                                if(measureMeanFreePath_)
                                {
                                    MFP_.boundaryField()[j][k] = 
                                                        MFP_[bCs[k]];
                                    SOF_.boundaryField()[j][k] =
                                                        SOF_[bCs[k]];
                                    MFPdX_.boundaryField()[j][k] =
                                                        MFPdX_[bCs[k]];
                                    MCR_.boundaryField()[j][k] =
                                                        MCR_[bCs[k]];
                                    MCT_.boundaryField()[j][k] =
                                                        MCT_[bCs[k]];
                                     MCTdt_.boundaryField()[j][k] =
                                                        MCTdt_[bCs[k]];
                                }
                                if(measureClassifications_)
                                {
                                    cIDist_.boundaryField()[j][k] = 
                                                        cIDist_[bCs[k]];
                                    cIIDist_.boundaryField()[j][k] = 
                                                        cIIDist_[bCs[k]];
                                    cIIIDist_.boundaryField()[j][k] = 
                                                        cIIIDist_[bCs[k]];
                                }
                                if(measureHeatFluxShearStress_)
                                {
                                    shearStressTensor_.boundaryField()[j][k] = 
                                                    shearStressTensor_[bCs[k]];
                                    heatFluxVector_.boundaryField()[j][k] = 
                                                        heatFluxVector_[bCs[k]];
                                    pressureTensor_.boundaryField()[j][k] = 
                                                        pressureTensor_[bCs[k]];
                                }

                                if(!isA<wallPolyPatch>(patch))
                                {
                                    translationalT_.boundaryField()[j][k] =
                                                    translationalT_[bCs[k]];
                                    rotationalT_.boundaryField()[j][k] =
                                                    rotationalT_[bCs[k]];
                                    vibrationalT_.boundaryField()[j][k] =
                                                    vibrationalT_[bCs[k]];
                                    overallT_.boundaryField()[j][k] = 
                                                    overallT_[bCs[k]];
                                    p_.boundaryField()[j][k] =
                                                            p_[bCs[k]];
                                    Ma_.boundaryField()[j][k] =
                                                            Ma_[bCs[k]];
                                    UMean_.boundaryField()[j][k] =
                                                            UMean_[bCs[k]];
                                }
                            }
                        }
                    }
                }
            }

            if(measureMeanFreePath_)
            {
                MFP_.write();
                MFPdX_.write();
                MCR_.write();
                MCT_.write();
                MCTdt_.write();
                MCS_.write();
                SOF_.write();
            }

            if(measureClassifications_)
            {
                cIDist_.write();
                cIIDist_.write();
                cIIIDist_.write();
            }

            if(measureErrors_)
            {
                densityError_.write();
                velocityError_.write();
                temperatureError_.write();
                pressureError_.write();
            }

            if(measureHeatFluxShearStress_)
            {
                heatFluxVector_.write();
                pressureTensor_.write();
                shearStressTensor_.write();
            }

            p_.write();
            translationalT_.write();
            rotationalT_.write();
            vibrationalT_.write();
            electronicT_.write();
            overallT_.write();
            q_.write();
            tau_.write();
            Ma_.write();
            UMean_.write();
            fD_.write();
        }

        //- reset
        if(time_.resetFieldsAtOutput())
        {
            nTimeSteps_ = 0.0;

            forAll(rhoNMean_, c)
            {
                rhoNMean_[c] = scalar(0.0);
                rhoMMean_[c] = scalar(0.0);
                linearKEMean_[c] = scalar(0.0);
                momentumMean_[c] = vector::zero;
                rotationalEMean_[c] = scalar(0.0);
                rotationalDofMean_[c] = scalar(0.0);
                rhoNMeanInt_[c] = scalar(0.0);
                molsElec_[c] = scalar(0.0),
                nClassI_[c] = scalar(0.0);
                nClassII_[c] = scalar(0.0);
                nClassIII_[c] = scalar(0.0);
                collisionSeparation_[c] = scalar(0.0);
                nColls_[c] = scalar(0.0);
                muu_[c] = scalar(0.0);
                muv_[c] = scalar(0.0);
                muw_[c] = scalar(0.0);
                mvv_[c] = scalar(0.0);
                mvw_[c] = scalar(0.0);
                mww_[c] = scalar(0.0);
                mcc_[c] = scalar(0.0);
                mccu_[c] = scalar(0.0);
                mccv_[c] = scalar(0.0);
                mccw_[c] = scalar(0.0);
                eu_[c] = scalar(0.0);
                ev_[c] = scalar(0.0);
                ew_[c] = scalar(0.0);
                e_[c] = scalar(0.0);
                rhoNMeanXnParticle_[c] = scalar(0.0);
                rhoMMeanXnParticle_[c] = scalar(0.0);
                momentumMeanXnParticle_[c] = vector::zero;
                linearKEMeanXnParticle_[c] = scalar(0.0);
            }

            forAll(electronicETotal_, iD)
            {
                forAll(electronicETotal_[iD], cell)
                {
                    electronicETotal_[iD][cell] = 0.0;
                    mccSpecies_[iD][cell] = 0.0;
                    nParcels_[iD][cell] = 0.0;
                    nGroundElectronicLevel_[iD][cell] = 0.0;
                    nFirstElectronicLevel_[iD][cell] = 0.0;
                    nParcelsXnParticle_[iD][cell] = 0.0;

                    forAll(vibrationalETotal_[iD], v)
                    {
                       vibrationalETotal_[iD][v][cell] = 0.0; 
                    }
                }
            }

            // reset boundary information
            forAll(rhoNBF_, j)
            {
                rhoNBF_[j] = 0.0;
                rhoMBF_[j] = 0.0;
                linearKEBF_[j] = 0.0;
                speciesRhoNIntBF_[j] = 0.0;
                speciesRhoNElecBF_[j] = 0.0;
                rotationalEBF_[j] = 0.0;
                rotationalDofBF_[j] = 0.0;
                qBF_[j] = 0.0;
                fDBF_[j] = vector::zero;
                momentumBF_[j] = vector::zero;
            }

            forAll(speciesRhoNBF_, i)
            {
                forAll(speciesRhoNBF_[i], j)
                { 
                    speciesRhoNBF_[i][j] = 0.0;
                    vibrationalEBF_[i][j] = 0.0;
                    electronicEBF_[i][j] = 0.0;
                    mccSpeciesBF_[i][j] = 0.0;
                }
            }
        }

        if(averagingAcrossManyRuns_)
        {
            writeOut();
        }
    }
}


//- write field
void dsmcVolFields::writeField()
{}

void dsmcVolFields::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    if (propsDict_.found("measureClassifications"))
    {
        measureClassifications_ = 
            Switch(propsDict_.lookup("measureClassifications"));

        if(measureClassifications_)
        {
            Info << "measureClassifications initiated." << endl;
        }
    }

    if (propsDict_.found("measureErrors"))
    {
        measureErrors_ = Switch(propsDict_.lookup("measureErrors"));

        if(measureErrors_)
        {
            Info << "measureErrors initiated." << endl;
        }
    }

    if (propsDict_.found("densityOnly"))
    {
        densityOnly_ = Switch(propsDict_.lookup("densityOnly"));

        if(densityOnly_)
        {
            Info << nl << "densityOnly initiated." << nl << endl;
        }
    }

    if (propsDict_.found("measureHeatFluxShearStress"))
    {
        measureHeatFluxShearStress_ = 
            Switch(propsDict_.lookup("measureHeatFluxShearStress"));

        if(measureHeatFluxShearStress_)
        {
            Info << "measureHeatFluxShearStress initiated." << endl;
        }
    }

    if(propsDict_.found("measureMeanFreePath"))
    {
        measureMeanFreePath_ = Switch(propsDict_.lookup("measureMeanFreePath"));

        if(measureMeanFreePath_)
        {
            Info << "measureMeanFreePath initiated." << endl;

            mfpReferenceTemperature_ = 
            readScalar(propsDict_.lookup("mfpReferenceTemperature"));
        }
    }

    if (propsDict_.found("averagingAcrossManyRuns"))
    {
        averagingAcrossManyRuns_ = 
            Switch(propsDict_.lookup("averagingAcrossManyRuns"));

        if(averagingAcrossManyRuns_)
        {
            Info << "averagingAcrossManyRuns initiated." << endl;
        }
    }

}


} // End namespace Foam

// ************************************************************************* //

