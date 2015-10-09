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
constantMoleculeProperties::constantMoleculeProperties
(
    const polyMesh& mesh,
    const reducedUnits& rU
)
:
    mesh_(mesh),
    rU_(rU)
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
    chargeNames_.setSize(N_);
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
    DynamicList<word> chargeSiteIdList; 
    
    // reading in information from constant/moleculeProperties dictionary 
    forAll(molList, i)
    {
        const entry& molEntryI = molList[i];
        const dictionary& subDict = molEntryI.dict();
        
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

            if(rU_.runReducedUnits())
            {
                siteRefPositions_[i][j] /= rU_.refLength();
            }
        }

        List<scalar> siteMasses = subDict.lookup("siteMasses");
        siteMasses_[i].setSize(siteMasses.size());
        
        forAll(siteMasses, j)
        {
            siteMasses_[i][j] = siteMasses[j]; 
            
            if(rU_.runReducedUnits())
            {
                siteMasses_[i][j] /= rU_.refMass();
            }
        }        

        List<scalar> siteCharges = subDict.lookup("siteCharges");
        siteCharges_[i].setSize(siteCharges.size());
        
        DynamicList<word> chargeNames;
        
        forAll(siteCharges, j)
        {
            if(siteCharges[j] != 0.0)
            {
                if(findIndex(chargeSiteIdList, siteIdNames[j]) == -1)
                {
                    chargeSiteIdList.append(siteIdNames[j]);
                }
                
                chargeNames.append(siteNames_[i][j]);
            }
            
            siteCharges_[i][j] = siteCharges[j];
            
            if(rU_.runReducedUnits())
            {
                siteCharges_[i][j] /= rU_.refCharge();
            }
        }
        
        chargeNames_[i].transfer(chargeNames);
    }

    nPairPotSites_ = pairPotentialSiteIdList.size();
    nSites_=siteIdList.size();
    nChargeSites_=chargeSiteIdList.size();

    
    siteIdList_.transfer(siteIdList);    
    pairPotSiteIdList_.transfer(pairPotentialSiteIdList);    
    chargeSiteIdList_.transfer(chargeSiteIdList);
    
    
    pairPotNamesToFullSites_.setSize(N_);

    forAll(pairPotNamesToFullSites_, i)
    {
        pairPotNamesToFullSites_[i].setSize(pairPotNames_[i].size());
        
        forAll(pairPotNames_[i], j)
        {
            const word& name = pairPotNames_[i][j];
            label k = findIndex(pairPotSiteIdList_, name);
            pairPotNamesToFullSites_[i][j] = k;
        }
    }

    pairPotNamesToSites_.setSize(N_);

    forAll(pairPotNamesToSites_, i)
    {
        pairPotNamesToSites_[i].setSize(pairPotNames_[i].size());
        
        forAll(pairPotNames_[i], j)
        {
            const word& name = pairPotNames_[i][j];
            label k = findIndex(siteNames_[i], name);
            pairPotNamesToSites_[i][j] = k;
        }
    }

    
    chargePotNamesToFullSites_.setSize(N_);

    forAll(chargePotNamesToFullSites_, i)
    {
        chargePotNamesToFullSites_[i].setSize(chargeNames_[i].size());
        
        forAll(chargeNames_[i], j)
        {
            const word& name = chargeNames_[i][j];
            label k = findIndex(chargeSiteIdList_, name);
            chargePotNamesToFullSites_[i][j] = k;
        }
    }
    
    Info << "siteNames_ = " << siteNames_ << endl;
    Info << "chargeNames_ = " << chargeNames_ << endl;
    Info << "chargeSiteIdList_ = " << chargeSiteIdList_ << endl;
    Info << "chargePotNamesToFullSites_ = " << chargePotNamesToFullSites_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantMoleculeProperties::~constantMoleculeProperties()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const label& constantMoleculeProperties::nMolTypes() const
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

const label& constantMoleculeProperties::nSiteTypes() const
{
    return nSites_;
}

const label& constantMoleculeProperties::nPairPotTypes() const
{
    return nPairPotSites_;
}

const label& constantMoleculeProperties::nChargeTypes() const
{
    return nChargeSites_;
}


const List<word>& constantMoleculeProperties::siteIds() const
{
    return siteIdList_;
}

const List<word>& constantMoleculeProperties::pairPotSiteIdList() const
{
    return pairPotSiteIdList_;
}

const List<word>& constantMoleculeProperties::chargeSiteIdList() const
{
    return chargeSiteIdList_;
}

        
const List<List<label> >& constantMoleculeProperties::pairPotNamesToFullSites() const
{
    return pairPotNamesToFullSites_;
}

const List<List<label> >& constantMoleculeProperties::pairPotNamesToSites() const
{
    return pairPotNamesToSites_;
}

const List<List<label> >& constantMoleculeProperties::chargePotNamesToFullSites() const
{
    return chargePotNamesToFullSites_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
