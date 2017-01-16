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

#include "agentFieldProperties.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agentFieldProperties, 0);
addToRunTimeSelectionTable(agentMeasurement, agentFieldProperties, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



// Construct from components
agentFieldProperties::agentFieldProperties
(
    Time& t,
    const polyMesh& mesh,
    agentCloud& cloud,
    const dictionary& dict
)
:
    agentMeasurement(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    agentIds_(),
    fieldName_(propsDict_.lookup("fieldName")),
    boundaryCells_(),
    N_(mesh_.nCells(), 0.0),
    mass_(mesh_.nCells(), 0.0),
    mom_(mesh_.nCells(), vector::zero),
    vel_(mesh_.nCells(), vector::zero),
    kE_(mesh_.nCells(), 0.0),    

    NField_
    (
        IOobject
        (
            "N_"+ fieldName_,
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless, 0.0)
    ),
    rhoNField_
    (
        IOobject
        (
            "rhoN_"+ fieldName_,
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless/dimArea, 0.0)
    ),
    rhoMField_
    (
        IOobject
        (
            "rhoM_"+ fieldName_,
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimMass/dimArea, 0.0)
    ),
    UField_
    (
        IOobject
        (
            "U_"+ fieldName_,
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("0.0", dimVelocity, vector::zero)
    ),
    momField_
    (
        IOobject
        (
            "mom_"+ fieldName_,
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("0.0", dimMass*dimVelocity, vector::zero)
    ),  
    kEField_
    (
        IOobject
        (
            "kE_"+ fieldName_,
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimMass*dimVelocity*dimVelocity, 0.0)
    ),     
    nTimeSteps_(0.0),
    resetAtOutput_(true)  
{
   
    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();

    boundaryCells_.setSize(mesh.boundaryMesh().size());
    
    resetAtOutput_ = Switch(propsDict_.lookup("resetAtOutput"));
    
    boundedBox bbMesh( mesh_.bounds().min(), mesh_.bounds().max() );
    
       
    dZ_ = bbMesh.span() & vector(0, 0, 1);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentFieldProperties::~agentFieldProperties()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void agentFieldProperties::createField()
{
    NField_.write();
    rhoNField_.write();
    rhoMField_.write();
    UField_.write();
    momField_.write();
    kEField_.write();
}

void agentFieldProperties::calculateField()
{
    nTimeSteps_ += 1.0;
    
    {
        IDLList<agent>::iterator mol(cloud_.begin());
        
        for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
        {
            if(findIndex(agentIds_, mol().id()) != -1)
            {
                const label& cell = mol().cell();
                
                N_[cell] += 1.0;
                mass_[cell] += mol().mass();
                mom_[cell] += mol().v()*mol().mass();
                vel_[cell] += mol().v();
                kE_[cell] += 0.5*magSqr(mol().v());
            }
        }
    }
    
    if(time_.outputTime())
    {
        forAll(N_, cell)
        {
            if(N_[cell] > VSMALL)
            {
                scalar V = mesh_.cellVolumes()[cell];
                
                scalar totalArea_ = V/dZ_;
                
//                 if(cell < 1)
//                 {
//                     Info << "cell area = " << totalArea_ << endl;
//                 }
                
                NField_[cell] = N_[cell]/(nTimeSteps_); 
                rhoNField_[cell] = N_[cell]/(totalArea_*nTimeSteps_);
                rhoMField_[cell] = mass_[cell]/(totalArea_*nTimeSteps_);
                UField_[cell] = vel_[cell]/nTimeSteps_;
                kEField_[cell] = kE_[cell]/nTimeSteps_;
                momField_[cell] = mom_[cell]/nTimeSteps_;                
                
            }
            else
            {
                NField_[cell] = 0.0;
                rhoNField_[cell] = 0.0;
                rhoMField_[cell] = 0.0;
                UField_[cell] = vector::zero;
                momField_[cell] = vector::zero;
                kEField_[cell] = 0.0;
            }
        }
        
        forAll(boundaryCells_, j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            
            if(isA<polyPatch>(patch))
            {
                if(!isA<emptyPolyPatch>(patch))
                {
                    if(!isA<cyclicPolyPatch>(patch))
                    {
                        forAll(boundaryCells_[j], k)
                        {       
                            NField_.boundaryField()[j][k] = NField_[boundaryCells_[j][k]];
                            rhoNField_.boundaryField()[j][k] = rhoNField_[boundaryCells_[j][k]];
                            rhoMField_.boundaryField()[j][k] = rhoMField_[boundaryCells_[j][k]];
                            UField_.boundaryField()[j][k] = UField_[boundaryCells_[j][k]];
                            momField_.boundaryField()[j][k] = momField_[boundaryCells_[j][k]];
                            kEField_.boundaryField()[j][k] = kEField_[boundaryCells_[j][k]];
                        }
                    }
                }
            }
        }
        
        NField_.boundaryField() = NField_.boundaryField().boundaryInternalField();
        
        rhoNField_.correctBoundaryConditions();
        rhoMField_.correctBoundaryConditions();
        
        if(resetAtOutput_)
        {
            nTimeSteps_ = 0.0;
            
            forAll(N_, c)
            {
                N_[c] = scalar(0.0);
                mass_[c] = scalar(0.0);
                mom_[c] = vector::zero;
                vel_[c] = vector::zero;
                kE_[c] = scalar(0.0);
            }
        }
        
    }
}

void agentFieldProperties::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {


        }
    }    
}



void agentFieldProperties::measureDuringForceComputation
(
    agent* molI,
    agent* molJ
)
{}

void agentFieldProperties::measureDuringForceComputationSite
(
    agent* molI,
    agent* molJ,
    label sI,
    label sJ
)
{}

// const propertyField& agentFieldProperties::fields() const
// {
//     return fields_;
// }

} // End namespace Foam

// ************************************************************************* //
