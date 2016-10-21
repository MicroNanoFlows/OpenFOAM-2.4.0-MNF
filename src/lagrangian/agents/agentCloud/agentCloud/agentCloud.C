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

\*----------------------------------------------------------------------------*/

#include "agentCloud.H"
#include "agentConfigurations.H"
#include "fvMesh.H"
//#include "polyMolsToDelete.H"
//#include "polyMappingModels.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<agent>, 0);
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //





/*
// NEW //
void Foam::agentCloud::checkForOverlaps()
{
    if(p_.checkPotentialOverlaps())
    {
        const scalar& potLim = p_.potentialEnergyLimit();
        
        Info<< nl << "Removing high energy overlaps, limit = "
            << potLim
            << ", from removalOrder list = " << p_.removalOrder()
            << endl;

        label initialSize = this->size();

        if (Pstream::parRun())
        {
            reduce(initialSize, sumOp<label>());
        }

        buildCellOccupancy();

        prepareInteractions();

        DynamicList<agent*> molsToDelete;
        DynamicList<label> molsToDeleteTNs;

        agent* molI = NULL;
        agent* molJ = NULL;


        {
            // Real-Real interactions
            const labelListList& dil = iL_.dil();

            forAll(dil, d)
            {
                forAll(cellOccupancy_[d],cellIMols)
                {
                    molI = cellOccupancy_[d][cellIMols];
                    label idI = molI->id();
                    bool molIDeleted = false;
                    label tNI=molI->trackingNumber();

                    forAll(dil[d], interactingCells)
                    {
                        List<agent*> cellJ =
                        cellOccupancy_[dil[d][interactingCells]];

                        forAll(cellJ, cellJMols)
                        {
                            molJ = cellJ[cellJMols];
                            label tNJ = molJ->trackingNumber();

                            label molJDeleted = findIndex(molsToDeleteTNs, tNJ);

                            if(!molIDeleted && (molJDeleted == -1))
                            {
                                if(evaluatePotentialLimit(molI, molJ, potLim))
                                {
                                    label idJ = molJ->id();

                                    label removeIdI = findIndex(p_.removalOrder(), idI);
                                    label removeIdJ = findIndex(p_.removalOrder(), idJ);

                                    if(removeIdI < removeIdJ)
                                    {
                                        molsToDelete.append(molI);
                                        molsToDeleteTNs.append(tNI);
                                        molIDeleted = true;
                                    }
                                    else
                                    {
                                        molsToDelete.append(molJ);
                                        molsToDeleteTNs.append(tNJ);
                                    }
                                }
                            }
                        }
                    }

                    forAll(cellOccupancy_[d], cellIOtherMols)
                    {
                        molJ = cellOccupancy_[d][cellIOtherMols];
                        label tNJ = molJ->trackingNumber();

                        label molJDeleted = findIndex(molsToDeleteTNs, tNJ);

                        if(!molIDeleted && (molJDeleted == -1))
                        {
                            if (molJ > molI)
                            {
                                if(evaluatePotentialLimit(molI, molJ, potLim))
                                {
                                    label idJ = molJ->id();

                                    label removeIdI = findIndex(p_.removalOrder(), idI);
                                    label removeIdJ = findIndex(p_.removalOrder(), idJ);

                                    if(removeIdI < removeIdJ)
                                    {
                                        molsToDelete.append(molI);
                                        molsToDeleteTNs.append(tNI);
                                        molIDeleted = true;
                                    }
                                    else
                                    {
                                        molsToDelete.append(molJ);
                                        molsToDeleteTNs.append(tNJ);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        {
            // Real-Referred interactions
            forAll(iL_.refCellsParticles(), r)
            {
                const List<label>& realCells = iL_.refCells()[r].neighbouringCells();

                forAll(iL_.refCellsParticles()[r], i)
                {
                    molJ = iL_.refCellsParticles()[r][i];
                    label tNJ = molJ->trackingNumber();

                    label molJDeleted = findIndex(molsToDeleteTNs, tNJ);

                    if(molJDeleted == -1)
                    {
                        forAll(realCells, rC)
                        {
                            List<agent*> molsInCell = cellOccupancy_[realCells[rC]];

                            forAll(molsInCell, j)
                            {
                                molI = molsInCell[j];
                                label tNI = molI->trackingNumber();

                                label molIDeleted = findIndex(molsToDeleteTNs, tNI);

                                if(molIDeleted == -1)
                                {
                                    if(evaluatePotentialLimit(molI, molJ, potLim))
                                    {
                                        label idJ = molJ->id();
                                        label idI = molI->id();

                                        label removeIdI = findIndex(p_.removalOrder(), idI);
                                        label removeIdJ = findIndex(p_.removalOrder(), idJ);

                                        if(removeIdI < removeIdJ)
                                        {
                                            molsToDelete.append(molI);
                                            molsToDeleteTNs.append(tNI);
                                        }
                                        else if(removeIdI == removeIdJ)
                                        {
                                            if (molI->trackingNumber() > molJ->trackingNumber())
                                            {
                                                molsToDelete.append(molI);
                                                molsToDeleteTNs.append(tNI);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        label nMolsDeleted = 0;

        forAll (molsToDelete, mTD)
        {
            nMolsDeleted++;

            Pout << nl << " WARNING: Deleting molecule "
                <<  " proc no = " << Pstream::myProcNo()
                << ", position = " << molsToDelete[mTD]->position()
                << ", molecule = " << cP_.molIds()[molsToDelete[mTD]->id()]
                << endl;

            deleteParticle(*(molsToDelete[mTD]));

        }

        if (Pstream::parRun())
        {
            reduce(nMolsDeleted, sumOp<label>());
        }

        if(nMolsDeleted > 0)
        {
            Info << nl << " WARNING: Total number of molecules deleted = " << nMolsDeleted << endl;
        }
        else
        {
            Info << " NO OVERLAPPING MOLECULES" << endl;
        }

        molsToDelete.clear();
    }
    
    buildCellOccupancy();

    prepareInteractions();
}

*/


void Foam::agentCloud::checkMoleculesInMesh()
{
    Info << nl << "checking cell-molecule addressing" << endl;

    DynamicList<agent*> molsToDelete;

    label initialSize = this->size();

    iterator mol(this->begin());

    label noOfModifiedMols = 0;

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        
        label cell = -1;
        label tetFace = -1;
        label tetPt = -1;

        mesh_.findCellFacePt
        (
            mol().position(),
            cell,
            tetFace,
            tetPt
        );
        
        if(cell != -1)
        {
            if(mol().cell() != cell)
            {
                mol().cell() = cell;
                mol().tetFace() = tetFace;
                mol().tetPt() = tetPt;
                noOfModifiedMols++;
            }
        }
        else
        {
            agent* molI = &mol();
            molsToDelete.append(molI);
        }
    }

    if(noOfModifiedMols > 0)
    {
        Pout<< tab << " molecules that changed cell = " 
            << noOfModifiedMols
            << endl;
    }

    forAll (molsToDelete, mTD)
    {
        Pout << nl << " WARNING: Molecule Outside Mesh - Deleting molecule "
                    <<  " proc no = " << Pstream::myProcNo()
                    << ", position = " << molsToDelete[mTD]->position()
                    << ", molecule = " << cP_.agentIds()[molsToDelete[mTD]->id()]
                    << endl;

        deleteParticle(*(molsToDelete[mTD]));
    }

    label molsRemoved = initialSize - this->size();

    if (Pstream::parRun())
    {
        reduce(molsRemoved, sumOp<label>());
    }

    Info<< tab <<" molecules removed from outside mesh = " 
        << molsRemoved 
        << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


//- Use for running agentFOAM
Foam::agentCloud::agentCloud
(
    Time& t,
    const polyMesh& mesh,
    const agentProperties& cP
//     cachedRandomMD& rndGen
)
:
    Cloud<agent>(mesh, "agentCloud", false),
    mesh_(mesh),
    cP_(cP),
//     rndGen_(rndGen),
    rndGen_(clock::getTime()),
    int_(t, mesh_, *this),
    rU_(1),
    cellOccupancy_(mesh_.nCells()),
    fields_(t, mesh_, *this),
    controllers_(t, mesh, *this),
    agentTracking_(),
    cyclics_(t, mesh_, -1),
    f_(t, mesh_, *this),
    b_(t, mesh_, *this),
    s_(t, mesh_, *this)
{
    agent::readFields(*this);

    checkMoleculesInMesh();

    // read in tracking numbers
    updateTrackingNumbersAfterRead();

    //check and remove high energy overalps
//     checkForOverlaps();
    
    buildCellOccupancy();
    
    int_.integrator()->init();
    
    
    // TESTS
//     writeReferredCloud();
}



//- general constructor
Foam::agentCloud::agentCloud
(
    Time& t,
    const polyMesh& mesh,
    const agentProperties& cP,
//     cachedRandomMD& rndGen, 
    const word& option,
    const bool& clearFields
)
    :
    Cloud<agent>(mesh, "agentCloud", false),
    mesh_(mesh),
    cP_(cP),
//     rndGen_(rndGen),    
    rndGen_(clock::getTime()),
    int_(t, mesh_, *this), 
    rU_(1),
    cellOccupancy_(mesh_.nCells()),
    fields_(t, mesh_),
    controllers_(t, mesh),
    agentTracking_(),
    cyclics_(t, mesh_, -1),
    f_(t, mesh_, *this),
    b_(t, mesh_, *this),
    s_(t, mesh_, *this)
{
    agent::readFields(*this);

    label initialParticles = this->size();

    if (Pstream::parRun())
    {
        reduce(initialParticles, sumOp<label>());
    }
   
    if(clearFields)
    {
        Info << "clearing existing field of particles  " << endl;

        clear();

        initialParticles = 0;
    }
    else
    {
        updateTrackingNumbersAfterRead();
    }

    if((option == "initialise") && clearFields)
    {
        agentConfigurations conf(mesh, *this);
        conf.setInitialConfig();
        buildCellOccupancy();
    }
    else if((option == "initialise") && !clearFields)
    {
        checkMoleculesInMesh();
        buildCellOccupancy();
        agentConfigurations conf(mesh, *this);
        conf.setInitialConfig();
    }
//     else if(option == "delete")
//     {
//         checkMoleculesInMesh();
//         buildCellOccupancy();
//         prepareInteractions();
//         polyMolsToDelete molsDel(mesh_, *this);
//     }
//     else if(option == "mapping")
//     {
//         polyMappingModels molsToMap(mesh_, *this);
//         buildCellOccupancy();
//     }
//     else if(option == "quickMapping")
//     {
//         checkMoleculesInMesh();
//         buildCellOccupancy();
//     }
    else if(option == "NULL")
    {
        buildCellOccupancy();
    }
    else 
    {
        Info << "ERROR" << endl;
    }

    label finalParticles = this->size();
    
    if (Pstream::parRun())
    {
        reduce(finalParticles, sumOp<label>());
    }

    
    Info << nl << "Initial agents = " << initialParticles 
         << ", modified agents = " << finalParticles - initialParticles
         << ", total agents: " << finalParticles 
         << endl;
}

// * * * * * * * * * * * * * * * * Static Constructors  * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::agentCloud> Foam::agentCloud::New
(
    Time& t,
    const polyMesh& mesh,
    const agentProperties& cP
//     cachedRandomMD& rndGen
)
{
    return autoPtr<agentCloud>
    (
        new agentCloud(t, mesh, cP)
    );
}

Foam::autoPtr<Foam::agentCloud> Foam::agentCloud::New
(
    Time& t,
    const polyMesh& mesh,
    const agentProperties& cP, 
//     cachedRandomMD& rndGen,
    const word& option,
    const bool& clearFields
)
{
    return autoPtr<agentCloud>
    (
        new agentCloud(t, mesh, cP, option, clearFields)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void  Foam::agentCloud::createAgent
(
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const vector& v,
    const vector& d,
    const vector& f,
    const vector& specialPosition,
    const scalar& mass,
    const scalar& potentialEnergy, 
    const scalar& R,
    const scalar& frac,   
    const label special,
    const label id,
    const label trackingNumber
)
{
    addParticle
    (
        new agent
        (
            mesh_,
            cP_,
            position,
            cellI,
            tetFaceI,
            tetPtI,
            v,
            d,
            f,
            specialPosition,
            mass,
            potentialEnergy,
            R,
            frac,
            special,
            id,
            trackingNumber
        )
    );
}


// Evolve functions 

void Foam::agentCloud::evolve()
{
    int_.integrator()->evolve();
}

// move molecules (tracking)
void Foam::agentCloud::move()
{
    agent::trackingData td1(*this, 1);
    Cloud<agent>::move(td1, mesh_.time().deltaTValue());
}


// Member Operators


Foam::label Foam::agentCloud::nAgents() const
{
    label size = this->size();

    if (Pstream::parRun())
    {
        reduce(size, sumOp<label>());
    }

    return size;
}

void Foam::agentCloud::updateTrackingNumbersAfterRead()
{
    const_iterator mol(this->begin());
    
    label tN = 0;
    
    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        if(mol().trackingNumber() > tN)
        {
            tN = mol().trackingNumber();
        }
    } 
    
    //- parallel-processing
    if(Pstream::parRun())
    {

        //- sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << tN;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                label trackingNumberProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> trackingNumberProc;
                }

                if(trackingNumberProc > tN)
                {
                    tN = trackingNumberProc;
                }
            }
        }
    }    
    
    agentTracking_.trackingIndex() = tN+1;
}

Foam::label Foam::agentCloud::getTrackingNumber()
{
    return agentTracking_.getTrackingNumber();
}

// This function is not being used.It is a test to see if the max number of possible 
// labels have been exceeded. The maximum is usually so large that this is hardly ever possible.
// It may only be required for systems with an enormous turnover of molecules (adds and deletes).   
void Foam::agentCloud::resetTrackingNumbers()
{
    agentTracking_.resetTrackingNumbers();

    if(agentTracking_.resetTracking())
    {
        iterator mol(this->begin());

        for (mol = this->begin(); mol != this->end(); ++mol)
        {
            mol().trackingNumber() = getTrackingNumber();
        }
    }
}

void Foam::agentCloud::buildCellOccupancy()
{
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    iterator mol(this->begin());

    for
    (
        mol = this->begin();
        mol != this->end();
        ++mol
    )
    {
        cellOccupancy_[mol().cell()].append(&mol());
    }
    
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].shrink();
    }
}

void Foam::agentCloud::insertMolInCellOccupancy(agent* mol)
{
    cellOccupancy_[mol->cell()].append(mol);
}

void Foam::agentCloud::removeMolFromCellOccupancy
(
    agent* molI
)
{
    DynamicList<agent*> updatedMolsInCell(0);

    const label& cellI = molI->cell();

    {
        const List<agent*>& molsInCell = cellOccupancy_[cellI];
    
        forAll(molsInCell, m)
        {
            agent* molJ = molsInCell[m];
    
            if(molI != molJ)
            {
                updatedMolsInCell.append(molJ);
            }
        }
    }

    cellOccupancy_[cellI].clear();
    cellOccupancy_[cellI].transfer(updatedMolsInCell);
}


void Foam::agentCloud::removeMolFromCellOccupancy
(
    const label& cellMolId,
    const label& cell
)
{
    DynamicList<agent*> molsInCell(0);

    forAll(cellOccupancy_[cell], c)
    {
        if(c != cellMolId)
        {
            molsInCell.append(cellOccupancy_[cell][c]);
        }
    }

    cellOccupancy_[cell].clear();
    cellOccupancy_[cell].transfer(molsInCell);
}


//- used if you want to read a new field at every time-step from an input file
//- e.g. to be used in a utility that computes measurements
// Used by reconstructXmol utility - reconstructPar does not produce the XMOL
// files after parallel processing
void Foam::agentCloud::readNewField()
{
    label initialSize = this->size();

    clear();

    IOPosition<Cloud<agent> > ioP(*this);

    if (ioP.headerOk())
    {
        ioP.readData(*this, false);
        ioP.close();
    }
    else
    {
        // WARNING
        WarningIn("readNewField()")
            << "Cannot read particle positions file " << nl
            << "    " << ioP.objectPath() << nl
            << "    assuming the initial cloud contains 0 particles." << endl;        
    }    
    
    particle::readFields(*this);

    agent::readFields(*this);

    if (this->size() != initialSize)
    {
        Info << "Changed agentCloud size, from: " 
                << initialSize << ", to: " << this->size() << endl;
    }

//     setSiteSizesAndPositions();    
}

void Foam::agentCloud::writeXMOL()
{

    {
        fileName path = mesh_.time().timePath()/cloud::prefix/"agentCloud.xmol";
        
        OFstream os(path);

        os << this->size() << nl << "agentCloud site positions in angstroms" << nl;

        const_iterator mol(this->begin());

        for (mol = this->begin(); mol != this->end(); ++mol)
        {
            os << cP_.agentIds()[mol().id()]
                << ' ' << mol().position().x()
                << ' ' << mol().position().y()
                << ' ' << mol().position().z()
                << nl;
        }
    }
    
}


// ************************************************************************* //
