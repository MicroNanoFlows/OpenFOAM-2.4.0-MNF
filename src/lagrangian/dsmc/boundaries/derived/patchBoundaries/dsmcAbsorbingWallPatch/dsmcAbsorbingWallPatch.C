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

#include "dsmcAbsorbingWallPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcAbsorbingWallPatch, 0);

addToRunTimeSelectionTable
(dsmcPatchBoundary, dsmcAbsorbingWallPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcAbsorbingWallPatch::dsmcAbsorbingWallPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeIds_(),
    absorptionProbs_()
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcAbsorbingWallPatch::~dsmcAbsorbingWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcAbsorbingWallPatch::initialConfiguration()
{
    
}

void dsmcAbsorbingWallPatch::calculateProperties()
{

}

void dsmcAbsorbingWallPatch::controlParticle
(dsmcParcel& p, dsmcParcel::trackingData& td)
{
    measurePropertiesBeforeControl(p);
    Random& rndGen = cloud_.rndGen();
    const label& iD = findIndex(typeIds_, p.typeId());
    
    if(iD != -1) //- particle might be absorbed
    {
        scalar absorbtionProbability = absorptionProbs_[iD];
        
        if(absorbtionProbability > rndGen.scalar01()) //- absorbed
        {
            //- delete the particle
             td.keepParticle = false;
        }
        else //- diffuse reflection
        {
            diffuseReflection(p);
        }   
    }
    else //- otherwise, it is treated as a diffuse reflection
    {
        diffuseReflection(p);
    }
}

void dsmcAbsorbingWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}


void dsmcAbsorbingWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();

}

void dsmcAbsorbingWallPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("temperature"));
    
    //  read in the type ids

    const List<word> molecules (propsDict_.lookup("absorbedIds"));

    if(molecules.size() == 0)
    {
        
        FatalErrorIn("dsmcAbsorbingWallPatch::setProperties()")
            << "Cannot have zero typeIds being absorbed." << nl << "in: "
            << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

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

    //  set the type ids

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            
            FatalErrorIn("dsmcAbsorbingWallPatch::setProperties()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
    
    const dictionary& absorptionProbabilitiesDict
    (
        propsDict_.subDict("absorptionProbabilities")
    );
    
    absorptionProbs_.clear();

    absorptionProbs_.setSize(typeIds_.size(), 0.0);

    forAll(absorptionProbs_, i)
    {
        absorptionProbs_[i] = readScalar
        (
            absorptionProbabilitiesDict.lookup(moleculesReduced[i])
        );
    }
}

void dsmcAbsorbingWallPatch::diffuseReflection(dsmcParcel& p)
{
    vector& U = p.U();

    scalar& ERot = p.ERot();
    
    labelList& vibLevel = p.vibLevel();
    
    label& ELevel = p.ELevel();

    label typeId = p.typeId();

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

    const scalar& T = temperature_;

    scalar mass = cloud_.constProps(typeId).mass();

    scalar rotationalDof = 
                    cloud_.constProps(typeId).rotationalDegreesOfFreedom();
    
    scalar vibrationalDof = 
                    cloud_.constProps(typeId).vibrationalDegreesOfFreedom();


    U =
        sqrt(physicoChemical::k.value()*T/mass)
       *(
            rndGen.GaussNormal()*tw1
          + rndGen.GaussNormal()*tw2
          - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))*nw
        );

       
    ERot = cloud_.equipartitionRotationalEnergy(T, rotationalDof);

    
    vibLevel = 
        cloud_.equipartitionVibrationalEnergyLevel(T, vibrationalDof, typeId);
   
    
    ELevel = cloud_.equipartitionElectronicLevel
                    (
                        T,
                        cloud_.constProps(typeId).degeneracyList(),
                        cloud_.constProps(typeId).electronicEnergyList(),
                        typeId
                    );   
    
    U += velocity_;
  
    measurePropertiesAfterControl(p, 0.0);
}

} // End namespace Foam

// ************************************************************************* //
