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

#include "dsmcStickingDiffuseWallPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcStickingDiffuseWallPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcStickingDiffuseWallPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcStickingDiffuseWallPatch::dsmcStickingDiffuseWallPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    stickingProbability_(),
    nStrikes_(0),
    strikeRate_(0.0),
    nSteps_(0)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcStickingDiffuseWallPatch::~dsmcStickingDiffuseWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcStickingDiffuseWallPatch::initialConfiguration()
{
    
}

void dsmcStickingDiffuseWallPatch::calculateProperties()
{

}

void dsmcStickingDiffuseWallPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
{
    nStrikes_++;
    
//     measurePropertiesBeforeControl(p);
    
    Random& rndGen(cloud_.rndGen());

    vector& U = p.U();

    scalar& ERot = p.ERot();
    
    labelList& vibLevel = p.vibLevel();

    label typeId = p.typeId();
    
    label& stuckToWall = p.stuckToWall();
    
    scalarField& wallTemperature = p.wallTemperature();
    
    vectorField& wallVectors = p.wallVectors();

    vector nw = p.normal();
    nw /= mag(nw);

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    // Wall tangential velocity (flow direction)
    vector Ut = U - U_dot_nw*nw;

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
    
    if(stickingProbability_ < rndGen.scalar01())
    {
//         measurePropertiesBeforeControl(p);
        
        const scalar& T = temperature_;

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

        ERot = cloud_.equipartitionRotationalEnergy(T, rotationalDof);
        
        vibLevel = cloud_.equipartitionVibrationalEnergyLevel(T, vibrationalDof, typeId);
        
        U += velocity_;
        
//         measurePropertiesAfterControl(p);
        
        
    }
    else
    {
        //Particle becomes stuck to wall
        stuckToWall = 1;
        
        measurePropertiesBeforeControl(p);
        
        scalar preIE = 0.0;
        vector preIMom = vector::zero;
        
        scalar mass = cloud_.constProps(typeId).mass();
        
        preIE = 0.5*mass*(U & U) + ERot + cloud_.constProps(typeId).electronicEnergyList()[p.ELevel()];
        
        forAll(p.vibLevel(), i)
        {
           preIE +=  p.vibLevel()[i]*cloud_.constProps(typeId).thetaV()[i]*physicoChemical::k.value();
        }
        
        preIMom = mass*U;
        
        wallTemperature[3] = preIE;
        wallVectors[3] = preIMom;
        
        const scalar& T = temperature_;
        
        U = SMALL*
            sqrt(physicoChemical::k.value()*T/mass)
        *(
                rndGen.GaussNormal()*tw1
            + rndGen.GaussNormal()*tw2
            - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))*nw
            );

//         U = vector::zero;
        
        wallTemperature[0] = temperature_;
        
        label wppIndex = patchId_;
        const polyPatch& wpp = mesh_.boundaryMesh()[wppIndex];
        label wppLocalFace = wpp.whichFace(p.face());
        
        wallTemperature[1] = wppIndex;
        wallTemperature[2] = wppLocalFace;
        
        wallVectors[0] = tw1;
        wallVectors[1] = tw2;
        wallVectors[2] = nw;
//         Info << "Stuck" << endl;
    }
}

void dsmcStickingDiffuseWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{	
}


void dsmcStickingDiffuseWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();
    
    nSteps_++;
    
    strikeRate_ = (nStrikes_*cloud_.nParticle())/(patchSurfaceArea_*mesh_.time().deltaTValue()*1000);
    
    Info << "strikeRate_ = " << strikeRate_ << endl;
    
    nStrikes_ = 0;

}

void dsmcStickingDiffuseWallPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("temperature"));
    stickingProbability_ = readScalar(propsDict_.lookup("stickingProbability"));
//     residenceTime_ = readScalar(propsDict_.lookup("residenceTime"));
}

} // End namespace Foam

// ************************************************************************* //
