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

#include "meanFreePath.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(meanFreePath, 0);

addToRunTimeSelectionTable(polyField, meanFreePath, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //






// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
meanFreePath::meanFreePath
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
    fieldName_(propsDict_.lookup("fieldName")),
    molIds_()
{
    dCol_ = readScalar(propsDict_.lookup("collisionDistance"));
    
    dCol_ /= molCloud_.redUnits().refLength();
   // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();
    
    testInitialTimeScale_ = false;

    if (propsDict_.found("testInitialTimeScale"))
    {
        testInitialTimeScale_ = Switch(propsDict_.lookup("testInitialTimeScale"));
    }
    
    nameFile1_ = "MFP_"+fieldName_+"_collectorOfFreePaths.txt";
    nameFile2_ = "MFP_"+fieldName_+"_freePaths.txt";
    nameFile3_ = "MFP_"+fieldName_+"_nCollisions.txt";
    
    nameFile4_ = "MFP_"+fieldName_+"_startCollisionPositions_X.txt";
    nameFile5_ = "MFP_"+fieldName_+"_startCollisionPositions_Y.txt";
    nameFile6_ = "MFP_"+fieldName_+"_startCollisionPositions_Z.txt";
    
    nameFile7_ = "MFP_"+fieldName_+"_endCollisionPositions_X.txt";
    nameFile8_ = "MFP_"+fieldName_+"_endCollisionPositions_Y.txt";
    nameFile9_ = "MFP_"+fieldName_+"_endCollisionPositions_Z.txt";
        
    nameFile10_ = "MFP_"+fieldName_+"_startCollisionTimes.txt";
    nameFile11_ = "MFP_"+fieldName_+"_endCollisionTimes.txt";
    
    deltaT_ = time_.deltaT().value();    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meanFreePath::~meanFreePath()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void meanFreePath::createField()
{
    label largestTN = 0;
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        label tNI = mol().trackingNumber();
        
        if(tNI > largestTN)
        {
            largestTN = tNI;
        }
    }
    
    if (Pstream::parRun())
    {
        //- sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << largestTN;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                label largestTNProc;
    
                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> largestTNProc;
                }
    
                if(largestTNProc > largestTN)
                {
                    largestTN = largestTNProc;
                }
            }
        }
    }
    
    Info << "number of tracking numbers = " << largestTN << endl;
    
    listSize_ = largestTN+1;
    
    freePaths_.setSize(listSize_, 0.0);
    Rold_.setSize(listSize_, 0.0);
    collectorOfFreePaths_.setSize(listSize_);
    nCollisions_.setSize(listSize_, 0.0);
    
    startCollPosX_.setSize(listSize_);
    startCollPosY_.setSize(listSize_);
    startCollPosZ_.setSize(listSize_);
    
    endCollPosX_.setSize(listSize_);
    endCollPosY_.setSize(listSize_);
    endCollPosZ_.setSize(listSize_);    
    
    startCollTimes_.setSize(listSize_);
    endCollTimes_.setSize(listSize_);
    
    setRolds();
    
    nMols_ = molCloud_.nMols();
    
    readFromStorage();
}

void meanFreePath::calculateField()
{}

void meanFreePath::setRolds()
{
    Rold_ = 0.0;
    
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            label tNI = mol().trackingNumber();
            
            Rold_[tNI] = mol().R();   
        }
    }    
    
    // parallel processing
    
    forAll(Rold_, i)
    {
        if(Pstream::parRun())
        {
            reduce(Rold_[i], sumOp<scalar>());
        }
    }
}


void meanFreePath::afterForce()
{
    // measure free paths on individual procs
    {
        List<scalar> freePaths(listSize_, 0.0);
        
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            label tNI = mol().trackingNumber();
            
            if(mol().R() > dCol_)
            {
                freePaths[tNI] += mag(mol().v()*deltaT_);
            }
        }
        
        // parallel comms between procs
        forAll(freePaths, i)
        {
            if(Pstream::parRun())
            {
                reduce(freePaths[i], sumOp<scalar>());
            }
            
            // set free paths on every proc correctly
            freePaths_[i] += freePaths[i];
        }
    }
    
    // check for collisions & record collision position
    {
        List<scalar> collided(listSize_, false);        
        List<scalar> freePaths(listSize_, 0.0);
        List<vector> positions(listSize_, vector::zero);
        
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            label tNI = mol().trackingNumber();
            
            if( (mol().R() < dCol_) && (Rold_[tNI] > dCol_) ) // collision criteria
            {
                collided[tNI] = true;
                freePaths[tNI] = freePaths_[tNI];
                positions[tNI] = mol().position();
            }
        }

        // parallel comms 
        forAll(freePaths, i)
        {
            if(Pstream::parRun())
            {
                reduce(collided[i], sumOp<bool>());
                reduce(freePaths[i], sumOp<scalar>());
                reduce(positions[i], sumOp<vector>());
            }
        }
        
        forAll(freePaths, i)
        {
            if(collided[i])
            {
                collectorOfFreePaths_[i].append(freePaths[i]);
                freePaths_[i] = 0.0;
                nCollisions_[i] += 1.0;
                
                endCollPosX_[i].append(positions[i].x());
                endCollPosY_[i].append(positions[i].y());
                endCollPosZ_[i].append(positions[i].z());
                
                endCollTimes_[i].append(time_.timeOutputValue());
            }
        }
    }
    
    // check for **start** of new free path & initial position
    {
        List<scalar> startedFreePath(listSize_, false);
        List<vector> positions(listSize_, vector::zero);
        
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            label tNI = mol().trackingNumber();
            
            if( (mol().R() > dCol_) && (Rold_[tNI] < dCol_) ) // collision criteria
            {
                startedFreePath[tNI] = true;
                positions[tNI] = mol().position();
            }
        }

        // parallel comms 
        forAll(startedFreePath, i)
        {
            if(Pstream::parRun())
            {
                reduce(startedFreePath[i], sumOp<bool>());
                reduce(positions[i], sumOp<vector>());
            }
        }
        
        forAll(startedFreePath, i)
        {
            if(startedFreePath[i])
            {
                startCollPosX_[i].append(positions[i].x());
                startCollPosY_[i].append(positions[i].y());
                startCollPosZ_[i].append(positions[i].z());
                
                startCollTimes_[i].append(time_.timeOutputValue());
            }
        }
    }     
    
    setRolds();
    
    if(testInitialTimeScale_)
    {
        bool allMolsCollided = true;
        
        forAll(nCollisions_, i)
        {
            if(nCollisions_[i] < 1.0)
            {
                allMolsCollided = false;
            }
        }
        
        if(allMolsCollided)
        {
            Info << " ALL MOLECULES HAVE COLLIDED AT LEAST ONCE " << endl;
        }
    }
}

void meanFreePath::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        writeToStorage();
    }
}

void meanFreePath::writeToStorage()
{
    // note all output/input is in reduced units - carefully convert!
    
    fileName pathName(time_.path()/time_.timeName()/"uniform"/"poly");
    
    {
        OFstream file(pathName/nameFile1_);

        if(file.good())
        {
            file << collectorOfFreePaths_ << endl;
        }
        else
        {
            FatalErrorIn("void meanFreePath::writeToStorage()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    
    {
        OFstream file(pathName/nameFile2_);

        if(file.good())
        {
            file << freePaths_ << endl;
        }
        else
        {
            FatalErrorIn("void meanFreePath::writeToStorage()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    
    {
        OFstream file(pathName/nameFile3_);

        if(file.good())
        {
            file << nCollisions_ << endl;
        }
        else
        {
            FatalErrorIn("void meanFreePath::writeToStorage()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }

    {
        OFstream file(pathName/nameFile4_);

        if(file.good())
        {
            file << startCollPosX_ << endl;
        }
        else
        {
            FatalErrorIn("void meanFreePath::writeToStorage()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    {
        OFstream file(pathName/nameFile5_);

        if(file.good())
        {
            file << startCollPosY_ << endl;
        }
        else
        {
            FatalErrorIn("void meanFreePath::writeToStorage()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    {
        OFstream file(pathName/nameFile6_);

        if(file.good())
        {
            file << startCollPosZ_ << endl;
        }
        else
        {
            FatalErrorIn("void meanFreePath::writeToStorage()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    {
        OFstream file(pathName/nameFile7_);

        if(file.good())
        {
            file << endCollPosX_ << endl;
        }
        else
        {
            FatalErrorIn("void meanFreePath::writeToStorage()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    {
        OFstream file(pathName/nameFile8_);

        if(file.good())
        {
            file << endCollPosY_ << endl;
        }
        else
        {
            FatalErrorIn("void meanFreePath::writeToStorage()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    {
        OFstream file(pathName/nameFile9_);

        if(file.good())
        {
            file << endCollPosZ_ << endl;
        }
        else
        {
            FatalErrorIn("void meanFreePath::writeToStorage()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }    
    
    {
        OFstream file(pathName/nameFile10_);

        if(file.good())
        {
            file << startCollTimes_ << endl;
        }
        else
        {
            FatalErrorIn("void meanFreePath::writeToStorage()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    {
        OFstream file(pathName/nameFile11_);

        if(file.good())
        {
            file << endCollTimes_ << endl;
        }
        else
        {
            FatalErrorIn("void meanFreePath::writeToStorage()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }     
}

void meanFreePath::readFromStorage()
{
    Info << "reading from storage" << endl;
    
    fileName pathName(time_.path()/time_.timeName()/"uniform"/"poly");
    
    {
        IFstream file(pathName/nameFile1_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> collectorOfFreePaths_;
        }
    }
    {
        IFstream file(pathName/nameFile2_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> freePaths_;
        }
    }    
    {
        IFstream file(pathName/nameFile3_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> nCollisions_;
        }
    }
    {
        IFstream file(pathName/nameFile4_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> startCollPosX_;
        }
    }
    {
        IFstream file(pathName/nameFile5_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> startCollPosY_;
        }
    }
    {
        IFstream file(pathName/nameFile6_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> startCollPosZ_;
        }
    }
    {
        IFstream file(pathName/nameFile7_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> endCollPosX_;
        }
    }
    {
        IFstream file(pathName/nameFile8_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> endCollPosY_;
        }
    }
    {
        IFstream file(pathName/nameFile9_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> endCollPosZ_;
        }
    }    
    {
        IFstream file(pathName/nameFile10_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> startCollTimes_;
        }
    }
    {
        IFstream file(pathName/nameFile11_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> endCollTimes_;
        }
    } 
    
//     Info << "collectorOfFreePaths = " << collectorOfFreePaths_ << endl;
//     Info << "freePaths = " << freePaths_ << endl;
//     Info << "nCollisions = " << nCollisions_ << endl;        
    
}

void meanFreePath::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void meanFreePath::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}


const propertyField& meanFreePath::fields() const
{
    return fields_;
}



} // End namespace Foam

// ************************************************************************* //
