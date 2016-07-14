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

The wall temperature is a function of the heat flux from incident parcels.

The radiation equation Q = epsilon*sigma*T^4 is used to calculate the wall temperature.

This class gives a value of temperature for each face on the patch.

\*---------------------------------------------------------------------------*/

#include "dsmcFixedHeatFluxFieldPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcFixedHeatFluxFieldPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcFixedHeatFluxFieldPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcFixedHeatFluxFieldPatch::dsmcFixedHeatFluxFieldPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    wallTemperature_
    (
        IOobject
        (
            "wallTemperature",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    EcTot_(nFaces_, 0.0),
    EcTotInc_(nFaces_, 0.0),
    EcTotSum_(nFaces_, 0.0),
    EcTotIncSum_(nFaces_, 0.0),
//     heatFlux_(nFaces_, 0.0),
    newWallTemperature_(nFaces_, 0.0),
    desiredHeatFlux_(),
    relaxationFactor_(),
//     alpha_(0.0),
    stepCounter_(0),
    nSamples_(),
    resetFieldsAtOutput_(false)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();

    forAll(newWallTemperature_, f)
    {
        newWallTemperature_[f] = temperature_;
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcFixedHeatFluxFieldPatch::~dsmcFixedHeatFluxFieldPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcFixedHeatFluxFieldPatch::initialConfiguration()
{}


void dsmcFixedHeatFluxFieldPatch::calculateProperties()
{
    stepCounter_++;

    const scalar deltaT = mesh_.time().deltaTValue();
    
    forAll(EcTotSum_, f)
    {
        EcTotSum_[f] += EcTot_[f];
        EcTotIncSum_[f] += EcTotInc_[f];
        EcTot_[f] = 0.0;
        EcTotInc_[f] = 0.0;
    }
    
    if(stepCounter_ >= nSamples_/*time_.time().outputTime()*/)
    {
        
        forAll(EcTotSum_, f)
        {
//             EcTotSum_[f] += EcTot_[f];
//             EcTotIncSum_[f] += EcTotInc_[f];
//             EcTot_[f] = 0.0;
//             EcTotInc_[f] = 0.0;
            
            if(fabs(EcTotSum_[f]) > VSMALL) // zero face temperature not allowed!
            {
                const label& faceI = faces_[f];
                scalar fA = mag(mesh_.faceAreas()[faceI]);

                scalar heatFlux = (cloud_.nParticle()*EcTotSum_[f])/(deltaT*nSamples_*fA);
                
                Pout << "heatFlux = " << heatFlux << endl;
                
                scalar energyFluxInc = (cloud_.nParticle()*EcTotIncSum_[f])/(deltaT*nSamples_*fA);
                
                scalar oldWallTemperature = newWallTemperature_[f];
                
                newWallTemperature_[f] = oldWallTemperature
                    *(1.0 + relaxationFactor_*((heatFlux - desiredHeatFlux_)/(fabs(desiredHeatFlux_) + 100.0)));
                    
                if(newWallTemperature_[f] < VSMALL)
                {
                    newWallTemperature_[f] = temperature_;
                }
                
                Pout << "newWallTemperature_[f] = " << newWallTemperature_[f] << endl;
            }                                    
        }
        
        label wppIndex = patchId_;
        
        forAll(newWallTemperature_, f)
        {
            if(newWallTemperature_[f] > VSMALL)
            {
                wallTemperature_.boundaryField()[wppIndex][f] = newWallTemperature_[f];
            }
        }
        
//         //- reset
//         if(resetFieldsAtOutput_)
//         {
            stepCounter_ = 0.0;
            EcTotSum_ = 0.0;
            EcTotIncSum_ = 0.0;
//         }
    }
}

void dsmcFixedHeatFluxFieldPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
{
    measurePropertiesBeforeControl(p);

    vector& U = p.U();

    scalar& ERot = p.ERot();
    
    label& vibLevel = p.vibLevel();

    label typeId = p.typeId();
    
    scalar m = cloud_.constProps(typeId).mass();

    scalar preIE = 0.5*m*(U & U) + ERot + vibLevel*physicoChemical::k.value()*cloud_.constProps(typeId).thetaV();
    
    label faceId = findIndex(faces_, p.face());

    vector nw = p.normal();
    nw /= mag(nw);

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    // Wall tangential velocity (flow direction)
    vector Ut = U - U_dot_nw*nw;

    Random& rndGen(cloud_.rndGen());

    while (mag(Ut) < SMALL)
    {
        // If the incident velocity is parallel to the face normal, no
        // tangential direction can be chosen.  Add a perturbation to the
        // incoming velocity and recalculate.

        U = vector
        (
            U.x()*(0.8 + 0.2*rndGen.scalar01()),
            U.y()*(0.8 + 0.2*rndGen.scalar01()),
            U.z()*(0.8 + 0.2*rndGen.scalar01())
        );

        U_dot_nw = U & nw;

        Ut = U - U_dot_nw*nw;
    }

    // Wall tangential unit vector
    vector tw1 = Ut/mag(Ut);

    // Other tangential unit vector
    vector tw2 = nw^tw1;

//     label faceId = findIndex(faces_, p.face());
    const scalar& T = newWallTemperature_[faceId];
    
//     Pout << "T = " << T << endl;

    scalar mass = cloud_.constProps(typeId).mass();

    scalar rotationalDof = cloud_.constProps(typeId).rotationalDegreesOfFreedom();
    
    scalar vibrationalDof = cloud_.constProps(typeId).vibrationalDegreesOfFreedom();

    U =
        sqrt(physicoChemical::k.value()*T/mass)
       *(
            rndGen.GaussNormal()*tw1
          + rndGen.GaussNormal()*tw2
          - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))*nw
        );

    U += velocity_;

    ERot = cloud_.equipartitionRotationalEnergy(T, rotationalDof);
    
    vibLevel = cloud_.equipartitionVibrationalEnergyLevel(T, vibrationalDof, typeId);

    measurePropertiesAfterControl(p, 0.0);
    
    scalar postIE = 0.5*m*(U & U) + ERot + vibLevel*physicoChemical::k.value()*cloud_.constProps(typeId).thetaV();
    
    EcTot_[faceId] += (preIE - postIE);
    
    EcTotInc_[faceId] += preIE;
}

void dsmcFixedHeatFluxFieldPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void dsmcFixedHeatFluxFieldPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

//     Info << "Updating properties in dsmcFixedHeatFluxFieldPatch" << endl;
    
    setProperties();

}

void dsmcFixedHeatFluxFieldPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("initialTemperature"));
    desiredHeatFlux_ = readScalar(propsDict_.lookup("desiredHeatFlux"));
    relaxationFactor_ = readScalar(propsDict_.lookup("relaxationFactor"));
    nSamples_ = readScalar(propsDict_.lookup("nSamples"));
    
//     Info << "resetAtOutput = " << resetFieldsAtOutput_ << endl;
    
//     resetFieldsAtOutput_ = propsDict_.lookup("resetAtOutput");
    
//     Info << "resetAtOutput = " << resetFieldsAtOutput_ << endl;

//     alpha_ = readScalar(propsDict_.lookup("energyAccommodationCoeff"));
// 
//     if(alpha_ < 0 || alpha_ > 1)
//     {
// 
//         FatalErrorIn("dsmcFixedHeatFluxFieldPatch::setProperties()")
//             << "The value of energyAccommodationCoeff should be between 0 and 1: " 
//             << alpha_ << nl 
//             << exit(FatalError);
//     }
}

} // End namespace Foam

// ************************************************************************* //
