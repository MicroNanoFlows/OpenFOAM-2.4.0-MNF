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

#include "agentBoundaries.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Null Constructor 
agentBoundaries::agentBoundaries
(
    Time& t,
    const polyMesh& mesh
)
:    
    time_(t),
    agentBoundariesDict_
    (
        IOobject
        (
            "boundariesDict",
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    nCyclicBoundaryModels_(0),
    nPatchBoundaryModels_(0),
    
    patchBoundaryList_(),
    patchBoundaryNames_(),
    patchBoundaryIds_(),
    pBFixedPathNames_(),
    patchBoundaryModels_(),
    patchToModelId_(mesh.boundaryMesh().size(), -1),
    
    
    cyclicBoundaryList_(),
    cyclicBoundaryNames_(),
    cyclicBoundaryIds_(),
    cMFixedPathNames_(),
    cyclicBoundaryModels_(),
    cyclicBoundaryToModelId_(mesh.boundaryMesh().size(), -1)
{}


//- Constructor
agentBoundaries::agentBoundaries
(
    Time& t,
    const polyMesh& mesh,
    agentCloud& cloud
)
:
    time_(t),
    agentBoundariesDict_
    (
        IOobject
        (
            "boundariesDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nCyclicBoundaryModels_(0),
    nPatchBoundaryModels_(0),
    
    patchBoundaryList_(agentBoundariesDict_.lookup("patchBoundaries")),
    patchBoundaryNames_(patchBoundaryList_.size()),
    patchBoundaryIds_(patchBoundaryList_.size()),
    pBFixedPathNames_(patchBoundaryList_.size()),
    patchBoundaryModels_(patchBoundaryList_.size()),
    patchToModelId_(mesh.boundaryMesh().size(), -1),    
    
    cyclicBoundaryList_(agentBoundariesDict_.lookup("cyclicBoundaries")),
    cyclicBoundaryNames_(cyclicBoundaryList_.size()),
    cyclicBoundaryIds_(cyclicBoundaryList_.size()),
    cMFixedPathNames_(cyclicBoundaryList_.size()),
    cyclicBoundaryModels_(cyclicBoundaryList_.size()),
    cyclicBoundaryToModelId_(mesh.boundaryMesh().size(), -1)
{
    Info << "Creating the boundary models: " << nl << endl;

    
    if( patchBoundaryModels_.size() > 0 )
    {
        forAll(patchBoundaryModels_, p)
        {
            const entry& boundaryI = patchBoundaryList_[p];
            const dictionary& boundaryIDict = boundaryI.dict();
    
            patchBoundaryModels_[p] = autoPtr<agentPatchBoundary>
            (
                agentPatchBoundary::New(t, mesh, cloud, boundaryIDict)
            );
    
            patchBoundaryNames_[p] = patchBoundaryModels_[p]->type();
            patchBoundaryIds_[p] = p;
            nPatchBoundaryModels_++;
        }
    }

    checkPatchBoundaryModels(mesh);    
    
    //- cyclic boundaries
    
    if( cyclicBoundaryModels_.size() > 0 )
    {
        forAll(cyclicBoundaryModels_, c)
        {
            const entry& boundaryI = cyclicBoundaryList_[c];
            const dictionary& boundaryIDict = boundaryI.dict();
    
            cyclicBoundaryModels_[c] = autoPtr<agentCyclicBoundary>
            (
                agentCyclicBoundary::New(t, mesh, cloud, boundaryIDict)
            );
    
            cyclicBoundaryNames_[c] = cyclicBoundaryModels_[c]->type();
            cyclicBoundaryIds_[c] = c;
            nCyclicBoundaryModels_++;
        }
    }    
    
    checkCyclicBoundaryModels(mesh);

    if(nPatchBoundaryModels_ > 0) 
    {
        // directory: case/boundaries
        fileName boundariesPath(time_.path()/"boundaries");

        if( !isDir(boundariesPath) )
        {
            mkDir(boundariesPath);
        }

        // directory: case/boundaries/poly
        fileName polyBoundariesPath(boundariesPath/"agent");

        if( !isDir(polyBoundariesPath) )
        {
            mkDir(polyBoundariesPath);
        }

        // directory: case/boundaries/poly/patchBoundaryModels
        fileName patchBoundaryModelsPath(polyBoundariesPath/"patchBoundaryModels");
    
        if (!isDir(patchBoundaryModelsPath))
        {
            mkDir(patchBoundaryModelsPath);    
        }

        forAll(patchBoundaryModels_, p)
        {
            if(patchBoundaryModels_[p]->writeInCase())
            {
                // directory: case/boundaries/poly/patchBoundaryModels/<patchBoundaryModel>
                fileName patchBoundaryModelPath(patchBoundaryModelsPath/patchBoundaryNames_[p]);

                if (!isDir(patchBoundaryModelPath))
                {
                    mkDir(patchBoundaryModelPath);    
                }
    
                const word& patchName = patchBoundaryModels_[p]->patchName();

                // directory: case/controllers/poly/patchBoundaryModels/<patchBoundaryModel>/<patchName>    
                fileName patchPath(patchBoundaryModelPath/patchName);
   
                if (!isDir(patchPath))
                {
                    mkDir(patchPath);    
                }
    
                pBFixedPathNames_[p] = patchPath;
            }
        }
    }    
    
    //- creating directories
    if(nCyclicBoundaryModels_ > 0) 
    {
        // directory: case/boundaries
        fileName boundariesPath(time_.path()/"boundaries");

        if( !isDir(boundariesPath) )
        {
            mkDir(boundariesPath);
        }

        // directory: case/boundaries/agent
        fileName agentBoundariesPath(boundariesPath/"agent");

        if( !isDir(agentBoundariesPath) )
        {
            mkDir(agentBoundariesPath);
        }

        // directory: case/boundaries/agent/cyclicBoundaryModels
        fileName cyclicBoundaryModelsPath(agentBoundariesPath/"cyclicBoundaryModels");
    
        if (!isDir(cyclicBoundaryModelsPath))
        {
            mkDir(cyclicBoundaryModelsPath);    
        }

        forAll(cyclicBoundaryModels_, c)
        {
            if(cyclicBoundaryModels_[c]->writeInCase())
            {
                // directory: case/boundaries/agent/cyclicBoundaryModels/<cyclicBoundaryModel>
                fileName cyclicBoundaryModelPath(cyclicBoundaryModelsPath/cyclicBoundaryNames_[c]);

                if (!isDir(cyclicBoundaryModelPath))
                {
                    mkDir(cyclicBoundaryModelPath);    
                }
                
                const word& patchName = cyclicBoundaryModels_[c]->patchName();    

                // directory: case/controllers/agent/cyclicBoundaryModels/<cyclicBoundaryModel>/<patchName>      
                fileName patchPath(cyclicBoundaryModelPath/patchName);
   
                if (!isDir(patchPath))
                {
                    mkDir(patchPath);    
                }
    
                cMFixedPathNames_[c] = patchPath;
            }
        }
    }
}

agentBoundaries::~agentBoundaries()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void agentBoundaries::checkCyclicBoundaryModels(const polyMesh& mesh)
{
    label nPolyPatches = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];
    
        if(isA<cyclicPolyPatch>(patch))
        {
            label patchIndex = patch.index();

            forAll(cyclicBoundaryModels_, c)
            {
                const label& patchId = cyclicBoundaryModels_[c]->patchId();
 
                if(patchIndex == patchId)
                {
                    nPolyPatches++;
                    cyclicBoundaryToModelId_[patchi] = c;
                }
            }
        }
    }

    if(nPolyPatches != nCyclicBoundaryModels_)
    {
        FatalErrorIn("agentBoundaries::checkBoundaryModels(const polyMesh& mesh)")
            << nl
            << " Number of cyclic boundary models = "  << nCyclicBoundaryModels_ 
            << " chosen in the boundaryiesDict are inconsistent." 
            << abort(FatalError);
    }
}


void agentBoundaries::checkPatchBoundaryModels(const polyMesh& mesh)
{
    //- check that all poly-patches defined within blockMeshDict,
    //  each have one model.

    label nPolyPatches = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];
    
        if
        (
            isA<polyPatch>(patch) &&
            !isA<cyclicPolyPatch>(patch) &&
            !isA<processorPolyPatch>(patch) &&
            !isA<wallPolyPatch>(patch) &&
            !isA<emptyPolyPatch>(patch) &&
            !isA<symmetryPolyPatch>(patch)
        )
        {
            nPolyPatches++;

            label patchIndex = patch.index();

            label nPatches = 0;

            forAll(patchBoundaryModels_, p)
            {
                const label& patchId = patchBoundaryModels_[p]->patchId();
 
                if(patchIndex == patchId)
                {
                    nPatches++;
                    patchToModelId_[patchi] = p;
                }
            }

            if(nPatches > 1)
            {
                FatalErrorIn("polyBoundaries::checkPatchBoundaryModels(const polyMesh& mesh)")
                    << nl
                    << " Only one patch boundary model per poly-patch, [name: "
                    << patch.name()
                    << "]. No of models chosen for this patch are: " 
                    << nPatches  << ", in " 
                    << mesh.time().system()/"polyBoundariesDict"
                    << abort(FatalError);
            }
        }
    }

//     Pout << "patchToModelId_: " << patchToModelId_ << endl;

    if(nPolyPatches != nPatchBoundaryModels_)
    {
        FatalErrorIn("polyBoundaries::checkPatchBoundaryModels(const polyMesh& mesh)")
            << nl
            << " Number of poly-patches = "  << nPolyPatches 
            << " in blockMeshDict, are not equal to the number of patch models = " 
            << nPatchBoundaryModels_  << ", defined in " 
            << mesh.time().system()/"polyBoundariesDict"
            << abort(FatalError);
    }
}



void agentBoundaries::initialConfig()
{
    forAll(patchBoundaryModels_, p)
    {
        patchBoundaryModels_[p]->initialConfiguration();
    }    
    
    forAll(cyclicBoundaryModels_, c)
    {
        cyclicBoundaryModels_[c]->initialConfiguration();
    }
}

void agentBoundaries::calculateProperties()
{
    forAll(patchBoundaryModels_, p)
    {
        patchBoundaryModels_[p]->calculateProperties();
    }

    forAll(cyclicBoundaryModels_, c)
    {
        cyclicBoundaryModels_[c]->calculateProperties();
    }

//     forAll(generalBoundaryModels_, g)
//     {
//         generalBoundaryModels_[g]->calculateProperties();
//     }
}
/*
void agentBoundaries::controlAfterMove()
{
//     forAll(cyclicBoundaryModels_, c)
//     {
//         cyclicBoundaryModels_[c]->controlAfterMove();
//     }
}*/

// impose model after calculation of forces
// void agentBoundaries::controlAfterForces()
// {
//     forAll(generalBoundaryModels_, g)
//     {
//         generalBoundaryModels_[g]->controlMols();
//     }
// }

//- output
void agentBoundaries::outputResults()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        fileName uniformTimePath(runTime.path()/runTime.timeName()/"uniform");
        
        forAll(patchBoundaryModels_, p)
        {
            patchBoundaryModels_[p]->output(pBFixedPathNames_[p], uniformTimePath);
        }
    }
}



// const label& agentBoundaries::nPatchBoundaryModels() const
// {
//     return nPatchBoundaryModels_;
// }

const label& agentBoundaries::nCyclicBoundaryModels() const
{
    return nCyclicBoundaryModels_;
}

// const label& agentBoundaries::nGeneralBoundaryModels() const
// {
//     return nGeneralBoundaryModels_;
// }


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
