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

#include "dsmcReflectiveParticleMembranePatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcReflectiveParticleMembranePatch, 0);

addToRunTimeSelectionTable(dsmcCyclicBoundary, dsmcReflectiveParticleMembranePatch, dictionary);


void dsmcReflectiveParticleMembranePatch::readProperties()
{
    p_ = (readScalar(propsDict_.lookup("reflectionProbability")));
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("temperature"));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcReflectiveParticleMembranePatch::dsmcReflectiveParticleMembranePatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcCyclicBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    p_(readScalar(propsDict_.lookup("reflectionProbability"))),
    velocity_(propsDict_.lookup("velocity")),
    temperature_(readScalar(propsDict_.lookup("temperature"))),
    nReflections_(0),
    nRejections_(0)

{
    writeInTimeDir_ = false;
    writeInCase_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcReflectiveParticleMembranePatch::~dsmcReflectiveParticleMembranePatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void dsmcReflectiveParticleMembranePatch::calculateProperties()
{
    label nReflections = nReflections_;
    label nRejections = nRejections_;

    if(Pstream::parRun())
    {
        reduce(nReflections, sumOp<label>());
        reduce(nRejections, sumOp<label>());
    }

    if(nRejections > 0)
    {
        Info<< "no Reflections: " << nReflections 
            << ", no Rejections: " << nRejections
            << " ratio relfections/(reflections+rejections):" 
            << scalar(nReflections)/scalar(nReflections+nRejections)
            << endl;
    }
}

void dsmcReflectiveParticleMembranePatch::initialConfiguration()
{}



void dsmcReflectiveParticleMembranePatch::controlMol
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    const label& faceI = p.face();

    vector nF = mesh_.faceAreas()[faceI];
    vector nw = p.normal();
    
    vector& U = p.U();
    
//     const vector& fC  = mesh_.faceCentres()[faceI];
//     const label f = findIndex(controlPatch(), faceI);

//     bool reflect = false;

//     label f = findIndex(faces_, faceI);
    label fA = findIndex(coupledFacesA_, faceI);
    label fB = findIndex(coupledFacesB_, faceI);
    
    nF /= mag(nF);
    nw /= mag(nw);
    
    scalar U_dot_nw = U & nw;
    vector Ut = U - U_dot_nw*nw;
    
    label typeId = p.typeId();
    
    const scalar& T = temperature_;

    scalar mass = cloud_.constProps(typeId).mass();

    Random& rndGen = cloud_.rndGen();

    scalar d = nF & U;

/*            Info<< "parcel to reflect at pos: " 
                << p.position() << ", nF: " << nF
                << " old velocity: " << p.U() 
                << " faceI: " << faceI
                << " fA: " << fA
                << " fB: " << fB
                << " fB: " << fB
                << endl; */   


    if(d > 0) // processor boundary
    {
        if(fA != -1)
        {
            scalar pRandom = rndGen.scalar01();

            if( pRandom <= p_ ) // reflect molecule
            {
                //scalar Un = U & nF;
    
                //U -= 2.0*Un*nF;
                
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

                U =
                    sqrt(physicoChemical::k.value()*T/mass)
                *(
                        rndGen.GaussNormal()*tw1
                    + rndGen.GaussNormal()*tw2
                    - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))*nw
                    );
                
                scalar Un = U & nF;
    
                U -= 2.0*Un*nF;
                
//                 p.position() -= Un*nF*0.01*mesh_.time().deltaTValue();
                
                U += velocity_;

                td.switchProcessor = false;

                nReflections_++;

//                 Pout<< "Reflected!!: p at pos: " 
//                     << p.position() << ", nF: " << nF
//                     << " tracking number: " << p.trackingNumber()
//                     << " new velocity: " << p.v()
//                     << endl;

            }
            else
            {
                nRejections_++;
            }
        }
    }
    else if (d < 0) // cyclic (non-processor boundary)
    {
        if(fB != -1)
        {    
//             Info<< " to reflect p at pos: " 
//                 << p.position() << ", nF: " << nF
//                 << " old velocity: " << p.U() 
//                 << " faceI: " << faceI
//                 << " fA: " << fA
//                 << " fB: " << fB
//                 << endl;

            scalar pRandom = rndGen.scalar01();

            if( pRandom <= p_ ) // reflect molecule
            {
                //scalar Un = U & nF;
    
                //U -= 2.0*Un*nF;
                
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

                U =
                    sqrt(physicoChemical::k.value()*T/mass)
                *(
                        rndGen.GaussNormal()*tw1
                    + rndGen.GaussNormal()*tw2
                    - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))*nw
                    );
                
                scalar Un = U & nF;
    
                U -= 2.0*Un*nF;
                
//                 p.position() -= Un*nF*0.01*mesh_.time().deltaTValue();
                
                U += velocity_;

                nReflections_++;

//                 Info<< "Reflected!!: p at pos: " 
//                     << p.position() 
//                     << " new velocity: " << p.U()
//                     << endl;
            }
            else
            {
                nRejections_++;
            }
        }
    }

}



void dsmcReflectiveParticleMembranePatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}


void dsmcReflectiveParticleMembranePatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    readProperties();
}




} // End namespace Foam

// ************************************************************************* //
