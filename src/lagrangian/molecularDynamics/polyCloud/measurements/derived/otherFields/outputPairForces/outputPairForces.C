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

#include "outputPairForces.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(outputPairForces, 0);

addToRunTimeSelectionTable(polyField, outputPairForces, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// void outputPairForces::setBoundBox
// (
//     const dictionary& propsDict,
//     boundedBox& bb,
//     const word& name 
// )
// {
//     const dictionary& dict(propsDict.subDict(name));
//     
//     vector startPoint = dict.lookup("startPoint");
//     vector endPoint = dict.lookup("endPoint");
//     bb.resetBoundedBox(startPoint, endPoint);
// }




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
outputPairForces::outputPairForces
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fields_(t, mesh, "dummy"),
    fieldName_(propsDict_.lookup("fieldName"))    
{
    measureInterForcesSites_ = true;
    
    // choose molecule ids to sample
    {
        // choose molecule ids to sample
        molIdsWall_.clear();
        
        selectIds ids
        (
            molCloud_.cP(),
            propsDict_,
            "molIdsWall"
            
        );

        molIdsWall_ = ids.molIds();
    }
    
    {
        // choose molecule ids to sample
        molIdsFluid_.clear();
        
        selectIds ids
        (
            molCloud_.cP(),
            propsDict_,
            "molIdsFluid"
            
        );

        molIdsFluid_ = ids.molIds();
    }
    
    nameFile1_ = "outputPairs_"+fieldName_+"_positions.txt";
    nameFile2_ = "outputPairs_"+fieldName_+"_energies.txt";    
    nameFile3_ = "outputPairs_"+fieldName_+"_forces.txt";    
    nameFile4_ = "outputPairs_"+fieldName_+"_nPairs.txt";    
    nameFile5_ = "outputPairs_"+fieldName_+"_referenceMols.txt";    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

outputPairForces::~outputPairForces()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void outputPairForces::createField()
{
    label N = molCloud_.moleculeTracking().getMaxTrackingNumber();
    
    Pout << "N = " << N << endl;
    
    positions_.setSize(N, vector::zero);
    energies_.setSize(N, 0.0);
    forces_.setSize(N, vector::zero);
    nPairs_.setSize(N, 0.0);
    fluidMols_.setSize(N, false);

    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIdsFluid_, mol().id()) != -1 )
        {
            label tN = mol().trackingNumber();
            fluidMols_[tN] = true;
        }
    }     
    
    if(Pstream::parRun())
    {
        //-sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << fluidMols_;
                }
            }
        }

        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {

                List<bool> fluidMolsProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> fluidMolsProc;
                }
                
                forAll(fluidMolsProc, i)
                {
                    if(fluidMolsProc[i])
                    {
                        fluidMols_[i] = true;
                    }
                }
            }
        }
    }
}




void outputPairForces::calculateField()
{

}

void outputPairForces::afterForce()
{}

void outputPairForces::writeField()
{
    if(time_.outputTime())
    {
        fileName timePath(time_.path()/time_.timeName()/"uniform");
        
        {    
            IDLList<polyMolecule>::iterator mol(molCloud_.begin());

            for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
            {
                if(findIndex(molIdsFluid_, mol().id()) != -1)
                {
                    label tN = mol().trackingNumber();
                    positions_[tN] = mol().position();
                }
            }
        }

        
        {
            OFstream file(timePath/nameFile1_);

            if(file.good())
            {
                file << positions_ << endl;
            }
            else
            {
                FatalErrorIn("void outputPairForces::writeField()")
                    << "Cannot open file " << file.name()
                    << abort(FatalError);
            }
        }
        {
            OFstream file(timePath/nameFile2_);

            if(file.good())
            {
                file << energies_ << endl;
            }
            else
            {
                FatalErrorIn("void outputPairForces::writeField()")
                    << "Cannot open file " << file.name()
                    << abort(FatalError);
            }
        }        
        
        {
            OFstream file(timePath/nameFile3_);

            if(file.good())
            {
                file << forces_ << endl;
            }
            else
            {
                FatalErrorIn("void outputPairForces::writeField()")
                    << "Cannot open file " << file.name()
                    << abort(FatalError);
            }
        }    
        {
            OFstream file(timePath/nameFile4_);

            if(file.good())
            {
                file << nPairs_ << endl;
            }
            else
            {
                FatalErrorIn("void outputPairForces::writeField()")
                    << "Cannot open file " << file.name()
                    << abort(FatalError);
            }
        }  
        
        {
            OFstream file(timePath/nameFile5_);

            if(file.good())
            {
                file << fluidMols_ << endl;
            }
            else
            {
                FatalErrorIn("void outputPairForces::writeField()")
                    << "Cannot open file " << file.name()
                    << abort(FatalError);
            }
        }
        
        //clear fields
        positions_ = vector::zero;
        nPairs_ = 0.0;
        energies_ = 0.0;
        forces_ = vector::zero;               
    }
    

}

void outputPairForces::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

    
void outputPairForces::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{
    if(time_.outputTime())
    {
        label idI = molI->id();
        label idJ = molJ->id();
        
        label idIF = findIndex(molIdsFluid_, molI->id());
        label idJF = findIndex(molIdsFluid_, molJ->id());

        label idIW = findIndex(molIdsWall_, molI->id());
        label idJW = findIndex(molIdsWall_, molJ->id());

        label tN = -1; // trackingNumber of water molecule
        
        if
        (
            ((idIW != -1) && (idJF != -1)) 
        )
        {
            tN = molJ->trackingNumber();
        }
        else if((idJW != -1) && (idIF != -1))
        {
            tN = molI->trackingNumber();        
        }    
        
        if(tN != -1)
        {
            label k = molCloud_.pot().pairPots().pairPotentialIndex(idI, idJ, sI, sJ);    
            vector rsIsJ = molI->sitePositions()[sI] - molJ->sitePositions()[sJ];
            scalar rsIsJMag = mag(rsIsJ);
            scalar pE = molCloud_.pot().pairPots().energy(k, rsIsJMag);
            vector force = (rsIsJ/rsIsJMag)*molCloud_.pot().pairPots().force(k, rsIsJMag);
            
            if(molI->referred() || molJ->referred())
            {
                nPairs_[tN] += 0.5;
                energies_[tN] += pE*0.5*0.5;

                if(molI->trackingNumber() == tN)
                {
                    forces_[tN] += force*0.5;
                }
                else
                {
                    forces_[tN] -= force*0.5;
                }
            }
            else
            {
                nPairs_[tN] += 1;
                energies_[tN] += pE*0.5;
                
                if(molI->trackingNumber() == tN)
                {
                    forces_[tN] += force;
                }
                else
                {
                    forces_[tN] -= force;
                }
            }
        }
    }
}


const propertyField& outputPairForces::fields() const
{
    return fields_;
}



} // End namespace Foam

// ************************************************************************* //
