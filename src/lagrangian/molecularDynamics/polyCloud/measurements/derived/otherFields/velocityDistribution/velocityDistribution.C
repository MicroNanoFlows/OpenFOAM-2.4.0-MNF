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

#include "velocityDistribution.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(velocityDistribution, 0);

addToRunTimeSelectionTable(polyField, velocityDistribution, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void velocityDistribution::setBoundBox
(
    const dictionary& propsDict,
    boundedBox& bb,
    const word& name 
)
{
    const dictionary& dict(propsDict.subDict(name));
    
    vector startPoint = dict.lookup("startPoint");
    vector endPoint = dict.lookup("endPoint");

    bb.resetBoundedBox(startPoint, endPoint);
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
velocityDistribution::velocityDistribution
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
   // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();
    
    nameFile1_ = "velDistr_"+fieldName_+"_raw.txt";
    
    setBoundBox
    (
        propsDict_,
        bb_,
        "samplingBox" 
    );
    
    scalar binWidth = readScalar(propsDict_.lookup("binWidth_m_per_s"));
    
//     binWidth /= molCloud_.redUnits().refVelocity();
    
    velDistr_(binWidth);
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

velocityDistribution::~velocityDistribution()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void velocityDistribution::createField()
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
    
//     Info << "number of tracking numbers = " << largestTN << endl;
    
    listSize_ = largestTN+1;
    
    readFromStorage();
}

void velocityDistribution::calculateField()
{
    {
        List<scalar> vel(listSize_, 0.0);
        
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(bb_.contains(mol().position()))
            {
                label tNI = mol().trackingNumber();
                vel[tNI] = mag(mol().v())*molCloud_.redUnits().refVelocity();
            }
        }
        
        // parallel comms between procs
        forAll(vel, i)
        {
            if(Pstream::parRun())
            {
                reduce(vel[i], sumOp<scalar>());
            }
           
            velDistr_.add(vel[i]);
        }
    }    
}



void velocityDistribution::afterForce()
{
    


}

void velocityDistribution::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        writeToStorage();
    }
}

void velocityDistribution::writeToStorage()
{
    // note all output/input is in reduced units - carefully convert!

    fileName pathName(time_.path()/time_.timeName()/"uniform"/"poly");
    

    {
        List< Pair<scalar> > histogram = velDistr_.raw();        
        
        OFstream file(pathName/nameFile1_);

        if(file.good())
        {
            file << histogram << endl;
        }
        else
        {
            FatalErrorIn("void velocityDistribution::writeToStorage()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }        
    }
    
    {
        List< Pair<scalar> > histogram = velDistr_.normalised();        

        OFstream file(pathName/"velDistr_"+fieldName_+".txt");

        if(file.good())
        {
            forAll(histogram, i)
            {
                file 
                    << histogram[i].first() << "\t"
                    << histogram[i].second() << "\t"
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void mfpZone::write()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
}

void velocityDistribution::readFromStorage()
{
//     Info << "reading from storage" << endl;
    
    fileName pathName(time_.path()/time_.timeName()/"uniform"/"poly");
    
    {
        List< Pair<scalar> > histogram;
        
        IFstream file(pathName/nameFile1_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> histogram;
            
            forAll(histogram, i)
            {
                if(histogram[i].second() > 0)
                {
//                     Info << "second = " << histogram[i].second() << endl;
                    
                    for (label j = 0; j < label(histogram[i].second()); j++)
                    {
                        velDistr_.add(histogram[i].first());
                    }
                }
            }
            
//             Info << "distribution check = " <<  velDistr_.raw() << endl;
        }
    }      
    
}

void velocityDistribution::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void velocityDistribution::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}


const propertyField& velocityDistribution::fields() const
{
    return fields_;
}



} // End namespace Foam

// ************************************************************************* //
