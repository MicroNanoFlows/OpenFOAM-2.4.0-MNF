/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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


\*---------------------------------------------------------------------------*/

#include "agentCloud.H"
#include "agent.H"
#include "Time.H"
#include "cyclicPolyPatch.H"
#include "meshTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::agent::move
(
    agent::trackingData& td,
    const scalar& trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;
    
    const polyMesh& mesh = td.cloud().pMesh();
    
    if (special_ != SPECIAL_FROZEN)
    {
        scalar tEnd = (1.0 - stepFraction())*trackTime;
        scalar dtMax = tEnd;

        vector Utracking = v_;        
        
        while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
        {
            // Apply correction to position for reduced-D cases
            meshTools::constrainToMeshCentre(mesh, position());

            Utracking = v_;

            // Apply correction to velocity to constrain tracking for
            // reduced-D cases
            meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);
        
            // set the lagrangian time-step
            scalar dt = min(dtMax, tEnd);

            dt *= trackToFace(position() + dt*Utracking, td, false);

            tEnd -= dt;
            stepFraction() = 1.0 - tEnd/trackTime;
        }
    }

    return td.keepParticle;
}

void Foam::agent::setAsReferred()
{
    special_ = 1;
}


void Foam::agent::transformProperties(const tensor& T)
{
    particle::transformProperties(T);
}


void Foam::agent::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);

    if (special_ == SPECIAL_TETHERED)
    {
        specialPosition_ += separation;
    }
}



bool Foam::agent::hitPatch
(
    const polyPatch&,
    trackingData&,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}

void Foam::agent::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    td.switchProcessor = true;
}

void Foam::agent::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices& tetIs
)
{
    // Use of the normal from tetIs is not required as
    // hasWallImpactDistance for a moleculeCloud is false.
    vector nw = normal();
    nw /= mag(nw);

    scalar vn = v_ & nw;

    // Specular reflection
    if (vn > 0)
    {
        v_ -= 2*vn*nw;
    }
}

void Foam::agent::hitPatch
(
    const polyPatch& pp,
    trackingData& td
)
{
    //-find which patch has been hit
    label patchIndex = pp.index();

    const label& patchModelId = td.cloud().ob().patchToModelIds()[patchIndex];

    // apply a boundary model when a molecule collides with this poly patch
    td.cloud().ob().patchBoundaryModels()[patchModelId]->controlMol(*this, td);
}

void Foam::agent::hitCyclicPatch
(
    const cyclicPolyPatch& cpp, 
    trackingData& td
)
{
    //-find which patch has been hit
    label patchIndex = cpp.index();
    
    const label& patchModelId = td.cloud().ob().cyclicBoundaryToModelIds()[patchIndex];
    
    td.cloud().ob().cyclicBoundaryModels()[patchModelId]->control(*this, cpp, td);
    
}


// ************************************************************************* //

