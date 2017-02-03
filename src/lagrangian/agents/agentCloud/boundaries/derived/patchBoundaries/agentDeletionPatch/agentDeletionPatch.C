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

#include "agentDeletionPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agentDeletionPatch, 0);

addToRunTimeSelectionTable(agentPatchBoundary, agentDeletionPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
agentDeletionPatch::agentDeletionPatch
(
    Time& t,
    const polyMesh& mesh,
    agentCloud& cloud,
    const dictionary& dict
)
:
    agentPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
//     elapsedTime_(0.0),
//     writeInterval_(readScalar(t.controlDict().lookup("writeInterval"))),
//     startTime_(t.startTime().value())
    N_(0)
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentDeletionPatch::~agentDeletionPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void agentDeletionPatch::initialConfiguration()
{}

void agentDeletionPatch::calculateProperties()
{}

void agentDeletionPatch::controlMol
(
    agent& mol,
    agent::trackingData& td
)
{
    if(findIndex(agentIds_, mol.id()) != -1) 
    {
//         scalar massI = mol().mass();
//         massFlux_ += massI;
//         molFlux_ += 1.0;
//         cumulMolFlux_ += 1.0;
//         cumulMassFlux_ += massI;
        N_++;
        
        td.keepParticle = false;

//         Pout << "delete particle" << endl;
    }
    else // reflect
    {
        const label& faceI = mol.face();
        vector nF = mesh_.faceAreas()[faceI];
        nF /= mag(nF);
        
        scalar Un = mol.v() & nF;

        if (Un > 0.0)
        {
            mol.v() -= 2.0*Un*nF;
        }
    }

//     Pout << "mass flux: " << massFlux_ << endl;
}

void agentDeletionPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    label N = N_;
    
    if (Pstream::parRun())
    {
        reduce(N, sumOp<label>());
    }
    
    scalarField timeField(1, t_.timeOutputValue());
    scalarField NField(1, N);

    writeTimeData
    (
        fixedPathName,
        "deletionPatch_"+patchName_+"_N.xy",
        timeField,
        NField,
        true
    );


}

// void agentDeletionPatch::updateProperties(const dictionary& newDict)
// {
//     //- the main properties should be updated first
//     updateBoundaryProperties(newDict);
// }



} // End namespace Foam

// ************************************************************************* //
