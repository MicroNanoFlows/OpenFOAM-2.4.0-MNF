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

\*---------------------------------------------------------------------------*/

#include "constantMoleculeProperties.H"
// #include "writeTimeData.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//-  Constructor
constantMoleculeProperties::constantMoleculeProperties(const polyMesh& mesh)
:
    mesh_(mesh)
{
    Info << nl << "Reading moleculeProperties dictionary." << endl;
    
    IOdictionary dict
    (
        IOobject
        (
            "moleculeProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    
    List<word> idList = dict.lookup("idList");    
    
    N_ = idList.size();
    
    names_.setSize(N_);
    cloudTypes_.setSize(N_);
    siteNames_.setSize(N_);
    pairPotNames_.setSize(N_);
    siteRefPositions_.setSize(N_);
    siteMasses_.setSize(N_);
    siteCharges_.setSize(N_);
    
    forAll(names_, i)
    {
        if(findIndex(names_, idList[i]) != -1)
        {
            FatalErrorIn("constantMoleculeProperties.C") << nl
                    << " You have defined more than one time the molecule id = "
                    << idList[i]
                    << " in moleculeProperties, idList = "
                    << idList
                    << ". All molecule ids need to be unique"
                    << nl << abort(FatalError);            
        }
        
        names_[i] = idList[i];        
    }

    Info << "Info: number of molecule type ids = " << N_ 
         << nl << "names = " << names_
         << endl;
         
    PtrList<entry> molList(dict.lookup("moleculeProperties"));

    if(molList.size() != N_)
    {
        FatalErrorIn("constantMoleculeProperties.C") << nl
                << " Number of molecules in idList = " << N_
                << ", does not match number of molecule entries"
                << " in moleculeProperties subdictionary = "
                << molList.size()
                << nl << abort(FatalError);
    }
    
    // building unique lists
    DynamicList<word> siteIdList;
    DynamicList<word> pairPotentialSiteIdList;    

    // reading in information from constant/moleculeProperties dictionary 
    forAll(molList, i)
    {
        const entry& molEntryI = molList[i];
        const dictionary& subDict = molEntryI.dict();
//         Info << "entry name keyword = " <<  molEntryI.keyword()<< endl;
        
        // WE ARE MISSING THE CHECK ON THE HEADING OF THE SUBDICTIONARY,
        // AND THAT IT MATCHES WITH: names_
        
        if (molEntryI.keyword() != names_[i])
        {
            FatalErrorIn("constantMoleculeProperties.C")
                << "The header of the subdictionary = "  << molEntryI.keyword()
                << ", with id number = " << i 
                << ", does not match with the molecule id type = " << names_[i]
                << " in the full idList = "
                << names_
                << "WARNING: The order of the contents in the idList needs to "
                << " match the headers of the moleculeProperties subdictionaries."
                << nl << abort(FatalError);
        }        
        
        const word cloudTypeI = subDict.lookup("cloudType");   
        cloudTypes_[i] = cloudTypeI;
        List<word> siteIdNames = subDict.lookup("siteIds");
        siteNames_[i].setSize(siteIdNames.size());
        
        forAll(siteIdNames, j)
        {
            if(findIndex(siteIdList, siteIdNames[j]) == -1)
            {
                siteIdList.append(siteIdNames[j]);
            }

            siteNames_[i][j] = siteIdNames[j];
            
        }
        
        List<word> pairIdNames = subDict.lookup("pairPotentialSiteIds");
        pairPotNames_[i].setSize(pairIdNames.size());
        
        forAll(pairIdNames, j)
        {
            if (findIndex(siteIdNames, pairIdNames[j]) == -1)
            {
                FatalErrorIn("constantMoleculeProperties.C")
                    << "site id = " << pairIdNames[j] 
                    << " in pairPotentialSiteIds is not in its original siteIds: "
                    << siteIdNames << nl << abort(FatalError);
            }
            
            if (findIndex(pairPotentialSiteIdList, pairIdNames[j]) == -1)
            {
                pairPotentialSiteIdList.append(pairIdNames[j]);
            }
            
            pairPotNames_[i][j] = pairIdNames[j];
        }

        List<vector> siteRefPositions = subDict.lookup("siteReferencePositions");
        siteRefPositions_[i].setSize(siteRefPositions.size());
        
        forAll(siteRefPositions, j)
        {
            siteRefPositions_[i][j] = siteRefPositions[j];
        }        

        List<scalar> siteMasses = subDict.lookup("siteMasses");
        siteMasses_[i].setSize(siteMasses.size());
        
        forAll(siteMasses, j)
        {
            siteMasses_[i][j] = siteMasses[j];
        }        

        List<scalar> siteCharges = subDict.lookup("siteCharges");
        siteCharges_[i].setSize(siteCharges.size());
        
        forAll(siteCharges, j)
        {
            siteCharges_[i][j] = siteCharges[j];
        }        
    }

    nPairPotIds_ = pairPotentialSiteIdList.size();

    // not sure what is happening here
    // why would there 
//     forAll(siteIdList, k)
//     {
//         const word& siteName = siteIdList[k];
// 
//         if (findIndex(pairPotentialSiteIdList, siteName) == -1)
//         {
//             pairPotentialSiteIdList.append(siteName);
//         }
//     }

//     siteIdList_.transfer(pairPotentialSiteIdList);    
    
    siteIdList_.transfer(siteIdList);    
    potSiteIdList_.transfer(pairPotentialSiteIdList);    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantMoleculeProperties::~constantMoleculeProperties()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const label& constantMoleculeProperties::nTypes() const
{
    return N_;
}

const List<word>& constantMoleculeProperties::molIds() const
{
    return names_;
}

const List<word>& constantMoleculeProperties::cloudTypes() const
{
    return cloudTypes_;
}

const List<word>& constantMoleculeProperties::siteNames(const word& idName) const
{
    label id = findIndex(names_, idName);
    
    if (id == -1)
    {
        FatalErrorIn("constantMoleculeProperties.C") << nl
                << " Something went wrong trying to access the molecule name = " 
                << idName
                << ", from the idList = " << names_
                << nl << abort(FatalError);        
    }

    return siteNames_[id];
}

const List<word>& constantMoleculeProperties::siteNames(const label& id) const
{
    return siteNames_[id];
}

const List<word>& constantMoleculeProperties::pairPotNames(const word& idName) const
{
    label id = findIndex(names_, idName);
    
    if (id == -1)
    {
        FatalErrorIn("constantMoleculeProperties.C") << nl
                << " Something went wrong trying to access the molecule name = " 
                << idName
                << ", from the idList = " << names_
                << nl << abort(FatalError);        
    }

    return pairPotNames_[id];
}

const List<word>& constantMoleculeProperties::pairPotNames(const label& id) const
{
    return pairPotNames_[id];
}


const List<vector>& constantMoleculeProperties::siteRefPositions(const word& idName) const
{
    label id = findIndex(names_, idName);
    
    if (id == -1)
    {
        FatalErrorIn("constantMoleculeProperties.C") << nl
                << " Something went wrong trying to access the molecule name = " 
                << idName
                << ", from the idList = " << names_
                << nl << abort(FatalError);        
    }

    return siteRefPositions_[id];
}

const List<vector>& constantMoleculeProperties::siteRefPositions(const label& id) const
{
    return siteRefPositions_[id];
}

const List<scalar>& constantMoleculeProperties::siteMasses(const word& idName) const
{
    label id = findIndex(names_, idName);
    
    if (id == -1)
    {
        FatalErrorIn("constantMoleculeProperties.C") << nl
                << " Something went wrong trying to access the molecule name = " 
                << idName
                << ", from the idList = " << names_
                << nl << abort(FatalError);        
    }

    return siteMasses_[id];
}

const List<scalar>& constantMoleculeProperties::siteMasses(const label& id) const
{
    return siteMasses_[id];
}

const List<scalar>& constantMoleculeProperties::siteCharges(const word& idName) const
{
    label id = findIndex(names_, idName);
    
    if (id == -1)
    {
        FatalErrorIn("constantMoleculeProperties.C") << nl
                << " Something went wrong trying to access the molecule name = " 
                << idName
                << ", from the idList = " << names_
                << nl << abort(FatalError);        
    }

    return siteCharges_[id];
}

const List<scalar>& constantMoleculeProperties::siteCharges(const label& id) const
{
    return siteCharges_[id];
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
