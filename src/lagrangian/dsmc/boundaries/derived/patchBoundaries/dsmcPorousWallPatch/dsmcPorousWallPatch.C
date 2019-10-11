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

#include "dsmcPorousWallPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcPorousWallPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcPorousWallPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcPorousWallPatch::dsmcPorousWallPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    diffuseFraction_(readScalar(propsDict_.lookup("diffuseFraction")))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcPorousWallPatch::~dsmcPorousWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcPorousWallPatch::initialConfiguration()
{
    
}

void dsmcPorousWallPatch::calculateProperties()
{

}

void dsmcPorousWallPatch::controlParticle(dsmcParcel& p, 
                                                dsmcParcel::trackingData& td)
{
    Random& rndGen(cloud_.rndGen());
    
    if(diffuseFraction_ > rndGen.scalar01())
    {
        measurePropertiesBeforeControl(p);
        
    //     scalar currentTime = cloud_.mesh().time().value();
        
    //     Info << "currentTime = " << currentTime << endl;

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

        
        vibLevel = cloud_.equipartitionVibrationalEnergyLevel(
                T, 
                vibrationalDof,
                typeId
                                                             );
    
        
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
    else
    {
        //const scalar deltaT = mesh_.time().deltaTValue();
        
        //Move very slightly into next cell
        //p.position() += p.U()*deltaT*1e-6;
       
        //Info << "1" << endl;
        
        label& cellI = p.cell();
        
        Info << "cellI Before = " << cellI << endl;
        Info << "p.face() = " << p.face() << endl;
        
        if (cellI == mesh_.faceOwner()[p.face()])
        {
//             cellI = mesh_.faceNeighbour()[p.face()];
            Info << "Owner" << endl;
            Info << "Neighbour = " << mesh_.faceNeighbour()[p.face()] << endl;
        }
        else if (cellI == mesh_.faceNeighbour()[p.face()])
        {
//             cellI = mesh_.faceOwner()[p.face()];
            Info << "Neighbour" << endl;
            Info << "Owner = " << mesh_.faceOwner()[p.face()] << endl;
        }
        else
        {
            FatalErrorIn("Particle::trackToFace(const vector&, TrackData&)")
                << "addressing failure" << abort(FatalError);
        }
        
        //Info << "cellI After = " << cellI << endl;
        
       // Info << "2" << endl;
    }
}

void dsmcPorousWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}


void dsmcPorousWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();

}

void dsmcPorousWallPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("temperature"));
}

} // End namespace Foam

// ************************************************************************* //
