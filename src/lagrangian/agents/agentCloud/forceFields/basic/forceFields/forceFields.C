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

#include "forceFields.H"
#include "agentCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forceFields::forceFields
(
    Time& t,
    const polyMesh& mesh,
    agentCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    time_(t),
    cloud_(cloud),
    cP_(cloud_.cP()),
    forceFieldsDict_
    (
        IOobject
        (
            "forceFieldsDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    rCut_(readScalar(forceFieldsDict_.lookup("rCut"))/ (cloud.redUnits().refLength()) ) ,
    nPairPotentials_(0),
    
    pairPotList_(forceFieldsDict_.lookup("pairs")),
    pairPotNames_(pairPotList_.size()),
    pairPotIds_(pairPotList_.size()),
    pairPotentials_(pairPotList_.size()),
    iL_(mesh, cloud.redUnits(), cloud.cyclics(), rCut_, "agent")
//     ipl_(mesh.nCells())        
{

    Info << nl << "Creating forceFields" << nl << endl;

    //- create forceFields

    if(pairPotentials_.size() > 0 )
    {
        forAll(pairPotentials_, i)
        {
            const entry& pairPotI = pairPotList_[i];
            const dictionary& forceFieldsIDict = pairPotI.dict();
            const word& headerName=pairPotI.keyword();
            
            Info << "headerName = " << headerName << endl;
            
            pairPotentials_[i] = autoPtr<pairPotential>
            (
                pairPotential::New(cloud_, headerName, forceFieldsIDict)
            );
    
            pairPotNames_[i] = pairPotentials_[i]->type();
            pairPotIds_[i] = i;
    
            nPairPotentials_++;
        }
        
        testPairPotentials();
    }
    
    //test
    
    // making directory
    pathName_ = mesh_.time().path()/"forceFields";

    if(isDir(pathName_))
    {
        rmDir(pathName_);
    }     
    
    mkDir(pathName_);     
}

forceFields::~forceFields()
{}

//- initial configuration
//- call this function after the agentCloud is completely initialised
void forceFields::initialConfig()
{
    
    forAll(pairPotentials_, i)
    {
        pairPotentials_[i]->initialConfiguration();        
        pairPotentials_[i]->writeTables(pathName_);
    }    
    
}

void forceFields::testPairPotentials()
{
    Info << "test for pair potentials" << endl;
    
    const List<word>& idList = cP_.agentIds();
    
    pairPotentialLookUp_.setSize(idList.size());
    
    forAll(pairPotentialLookUp_, i)
    {
        pairPotentialLookUp_[i].setSize(idList.size(), -1);
    }
    
    //checking that all combination have been chosen
    
    forAll(pairPotentials_, i)
    {
        const List<word>& pairIdList = pairPotentials_[i]->idList();
        
        label idA = findIndex(idList, pairIdList[0]);
        label idB = findIndex(idList, pairIdList[1]);
        
        if(pairPotentialLookUp_[idA][idB] == -1)
        {
            pairPotentialLookUp_[idA][idB] = i;
        }
        else
        {
            FatalError
                << "forceFields::testPairPotentials(): " << nl
                << "    ids = " << pairIdList
                << " are already defined in system/forceFieldsDict"
                << endl;
        }
        
        if(pairPotentialLookUp_[idB][idA] == -1)
        {
            pairPotentialLookUp_[idB][idA] = i;
        }
        else
        {
            FatalError
                << "forceFields::testPairPotentials(): " << nl
                << "    ids = " << pairIdList
                << " are already defined in system/forceFieldsDict"
                << endl;
        }
    }
    
    forAll(pairPotentialLookUp_, i)
    {
        forAll(pairPotentialLookUp_[i], j)
        {
            if(pairPotentialLookUp_[i][j] == -1)
            {
                FatalError
                    << "forceFields::testPairPotentials(): " << nl
                    << "    ids = " << idList[i] << " - " << idList[j]
                    << " have not been defined in system/forceFieldsDict"
                    << endl;                
            }
        }
    }
    
    Info << "done" << endl;
}

void forceFields::calculatePairForces()
{
    iL_.setReferredParticles(cloud_.cellOccupancy());
    
    agent* molI = NULL;
    agent* molJ = NULL;

    {
        // Real-Real interactions
        const labelListList& dil = iL_.dil();
        const List<DynamicList<agent*> >& cO = cloud_.cellOccupancy();
        
        forAll(dil, d)
        {
            forAll(cO[d],cellIMols)
            {
                molI = cO[d][cellIMols];

                forAll(dil[d], interactingCells)
                {
                    List<agent*> cellJ = cO[dil[d][interactingCells]];

                    forAll(cellJ, cellJMols)
                    {
                        molJ = cellJ[cellJMols];
                            
                        evaluatePair(molI, molJ);
                    }
                }

                forAll(cO[d], cellIOtherMols)
                {
                    molJ = cO[d][cellIOtherMols];

                    if (molJ > molI)
                    {
                        evaluatePair(molI, molJ);
                    }
                }
            }
        }

        // Real-Referred interactions
        forAll(iL_.refCellsParticles(), r)
        {
            const List<label>& realCells = iL_.refCells()[r].neighbouringCells();

            forAll(iL_.refCellsParticles()[r], i)
            {
                molJ = iL_.refCellsParticles()[r][i];

                forAll(realCells, rC)
                {
                    List<agent*> molsInCell = cO[realCells[rC]];

                    forAll(molsInCell, j)
                    {
                        molI = molsInCell[j];
                        evaluatePair(molI, molJ);
                    }
                }
            }
        }
    }
}


/*
void forceFields::controlDuringForceComputation
(
    agent* molI,
    agent* molJ
)
{
    forAll(controllersDuringForceComp_, n)
    {
        const label& sC = controllersDuringForceComp_[n];
        stateControllers_[sC]->controlDuringForces(molI, molJ);
    }
}
*/



//- output -- call this function at the end of the MD time-step
// void forceFields::outputStateResults() 
// {
//  
// }

// const label& forceFields::nStateControllers() const
// {
//     return nStateControllers_;
// }


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
