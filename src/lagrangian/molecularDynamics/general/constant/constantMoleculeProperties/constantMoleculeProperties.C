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
        
        label Nsites = siteIdNames.size();
        
        siteNames_[i].setSize(Nsites);
        
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
        
        if(siteRefPositions.size() != Nsites)
        {
            FatalErrorIn("constantMoleculeProperties.C") << nl
                    << " Number of site reference positions = " << siteRefPositions.size()
                    << ", does not match number of sites "
                    << " in moleculeProperties subdictionary = "
                    << Nsites
                    << nl << abort(FatalError);
        }
    
        siteRefPositions_[i].setSize(Nsites, vector::zero);
        
        forAll(siteRefPositions, j)
        {
            siteRefPositions_[i][j] = siteRefPositions[j];

            if(rU_.runReducedUnits())
            {
                siteRefPositions_[i][j] /= rU_.refLength();
            }
        }

        List<scalar> siteMasses = subDict.lookup("siteMasses");
        
        if(siteMasses.size() != Nsites)
        {
            FatalErrorIn("constantMoleculeProperties.C") << nl
                    << " Number of siteMasses = " << siteMasses.size()
                    << ", does not match number of sites "
                    << " in moleculeProperties subdictionary = "
                    << Nsites
                    << nl << abort(FatalError);
        }
        
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
        
        if(siteCharges.size() != Nsites)
        {
            FatalErrorIn("constantMoleculeProperties.C") << nl
                    << " Number of siteCharges = " << siteCharges.size()
                    << ", does not match number of sites "
                    << " in moleculeProperties subdictionary = "
                    << Nsites
                    << nl << abort(FatalError);
        }
        
        siteCharges_[i].setSize(siteCharges.size(), 0.0);
        
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
    
    nSites_=siteIdList.size();
    nPairPotSites_ = pairPotentialSiteIdList.size();
    nChargeSites_=chargeSiteIdList.size();

    
    siteIdList_.transfer(siteIdList);    
    pairPotSiteIdList_.transfer(pairPotentialSiteIdList);    
    chargeSiteIdList_.transfer(chargeSiteIdList);
    
    
    pairPotNamesToPairPotSitesList_.setSize(N_);

    forAll(pairPotNamesToPairPotSitesList_, i)
    {
        pairPotNamesToPairPotSitesList_[i].setSize(pairPotNames_[i].size());
        
        forAll(pairPotNames_[i], j)
        {
            const word& name = pairPotNames_[i][j];
            label k = findIndex(pairPotSiteIdList_, name);
            pairPotNamesToPairPotSitesList_[i][j] = k;
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
    
    chargePotNamesToChargePotSitesList_.setSize(N_);

    forAll(chargePotNamesToChargePotSitesList_, i)
    {
        chargePotNamesToChargePotSitesList_[i].setSize(chargeNames_[i].size());
        
        forAll(chargeNames_[i], j)
        {
            const word& name = chargeNames_[i][j];
            label k = findIndex(chargeSiteIdList_, name);
            chargePotNamesToChargePotSitesList_[i][j] = k;
        }
    }
    
    chargePotNamesToSites_.setSize(N_);

    forAll(chargePotNamesToSites_, i)
    {
        chargePotNamesToSites_[i].setSize(chargeNames_[i].size());
        
        forAll(chargeNames_[i], j)
        {
            const word& name = chargeNames_[i][j];
            label k = findIndex(siteNames_[i], name);
            chargePotNamesToSites_[i][j] = k;
        }
    }    

    // Output tests - debug
    Info << "siteNames_ = " << siteNames_ << endl;
    Info << "pairPotNames_ = " << pairPotNames_ << endl;    
    Info << "chargeNames_ = " << chargeNames_ << endl;
    Info << "siteIdList_ = " << siteIdList_ << endl;
    Info << "pairPotSiteIdList_" << pairPotSiteIdList_ << endl;
    Info << "chargeSiteIdList_ = " << chargeSiteIdList_ << endl;


    // set mass
    
    masses_.setSize(N_, 0.0);
    
    
    forAll(siteMasses_, i)
    {
        forAll(siteMasses_[i], j)
        {
            masses_[i] += siteMasses_[i][j];
        }
    }
    
// apply adjustments and calculate moment of intertia
    
    momentOfInertia_.setSize(N_);
    
    forAll(names_, i)
    {
        vector centreOfMass(vector::zero);

        // Calculate the centre of mass of the body and subtract it from each
        // position

        forAll(siteRefPositions_[i], s)
        {
            centreOfMass += siteRefPositions_[i][s]*siteMasses_[i][s];
        }

        centreOfMass /= masses_[i];

        forAll(siteRefPositions_[i], s)
        {
            siteRefPositions_[i][s] -= centreOfMass;
        }

        if (siteNames_[i].size() == 1)
        {
            siteRefPositions_[i] = vector::zero;

            momentOfInertia_[i] = diagTensor(-1, -1, -1);
        }
        else if (linearMoleculeTest(i))
        {
            // Linear molecule.

            Info<< nl << "Linear molecule." << endl;

            vector dir = siteRefPositions_[i][1] -  siteRefPositions_[i][0];

            dir /= mag(dir);

            tensor Q = rotationTensor(dir, vector(1,0,0));

            // Transform the site positions
            forAll(siteRefPositions_[i], s)
            {
                siteRefPositions_[i][s] = (Q & siteRefPositions_[i][s]);
            }

            // The rotation was around the centre of mass but remove any
            // components that have crept in due to floating point errors

            centreOfMass = vector::zero;

            forAll(siteRefPositions_[i], s)
            {
                centreOfMass += siteRefPositions_[i][s]*siteMasses_[i][s];
            }

            centreOfMass /= masses_[i];

            forAll(siteRefPositions_[i], s)
            {
                siteRefPositions_[i][s] -= centreOfMass;
            }

            diagTensor momOfInertia = diagTensor::zero;

            forAll(siteRefPositions_[i], s)
            {
                const vector& p(siteRefPositions_[i][s]);
        
                momOfInertia += siteMasses_[i][s]*diagTensor(0, p.x()*p.x(), p.x()*p.x());
            }

            momentOfInertia_[i] = diagTensor
            (
                -1,
                momOfInertia.yy(),
                momOfInertia.zz()
            );
        }
        else
        {
            // Fully 6DOF molecule

            // Calculate the inertia tensor in the current orientation

            tensor momOfInertia(tensor::zero);
        
            forAll(siteRefPositions_[i], s)
            {
                const vector& p(siteRefPositions_[i][s]);
        
                momOfInertia += siteMasses_[i][s]*tensor
                (
                    p.y()*p.y() + p.z()*p.z(), -p.x()*p.y(), -p.x()*p.z(),
                    -p.y()*p.x(), p.x()*p.x() + p.z()*p.z(), -p.y()*p.z(),
                    -p.z()*p.x(), -p.z()*p.y(), p.x()*p.x() + p.y()*p.y()
                );
            }
            
    //         Info << "Eigen values: x = " << eigenValues(momOfInertia).x()
    //             << " y = " << eigenValues(momOfInertia).y()
    //             << " z = " << eigenValues(momOfInertia).z()
    //             << endl;
            Info << "momOfInertia (full tensor) = " <<  momOfInertia << endl;
            
            if (eigenValues(momOfInertia).x() < VSMALL)
            {
    //             FatalErrorIn("polyMolecule::constantProperties::constantProperties")
    //                 << "An eigenvalue of the inertia tensor is zero, the molecule "
    //                 << dict.name() << " is not a valid 6DOF shape."
    //                 << nl << abort(FatalError);
                
                Info << "Warning no adjustment to be made to molecule" << endl;
                
                momentOfInertia_[i] = diagTensor
                (
                    momOfInertia.xx(),
                    momOfInertia.yy(),
                    momOfInertia.zz()
                ); 
            }
            else
            {
                Info << "Adjusting molecule" << endl;
                
                // Normalise the inertia tensor magnitude to avoid SMALL numbers in the
                // components causing problems
        
                momOfInertia /= eigenValues(momOfInertia).x();
            
                tensor e = eigenVectors(momOfInertia);
        
            
                // Calculate the transformation between the principle axes to the
                // global axes
            
                tensor Q = vector(1,0,0)*e.x() + vector(0,1,0)*e.y() + vector(0,0,1)*e.z();
            
            //         Info << "Q: " << Q << endl;
            
                // Transform the site positions
                forAll(siteRefPositions_[i], s)
                {
                    siteRefPositions_[i][s] = (Q & siteRefPositions_[i][s]);
                }

                // Recalculating the moment of inertia with the new site positions
            
                // The rotation was around the centre of mass but remove any
                // components that have crept in due to floating point errors
            
                centreOfMass = vector::zero;
            
                forAll(siteRefPositions_[i], s)
                {
                    centreOfMass += siteRefPositions_[i][s]*siteMasses_[i][s];
                }
            
                centreOfMass /= masses_[i];
            
                forAll(siteRefPositions_[i], s)
                {
                    siteRefPositions_[i][s] -= centreOfMass;
            
    /*                Info<< "mol: " << sites_[s].name() 
                        << " position: " << sites_[s].siteReferencePosition()
                        << endl;*/
                }
                
                // Calculate the moment of inertia in the principle component
                // reference frame
            
                momOfInertia = tensor::zero;
            
                forAll(siteRefPositions_[i], s)
                {
                    const vector& p(siteRefPositions_[i][s]);
            
                    momOfInertia += siteMasses_[i][s]*tensor
                    (
                        p.y()*p.y() + p.z()*p.z(), -p.x()*p.y(), -p.x()*p.z(),
                        -p.y()*p.x(), p.x()*p.x() + p.z()*p.z(), -p.y()*p.z(),
                        -p.z()*p.x(), -p.z()*p.y(), p.x()*p.x() + p.y()*p.y()
                    );
                }
            
                momentOfInertia_[i] = diagTensor
                (
                    momOfInertia.xx(),
                    momOfInertia.yy(),
                    momOfInertia.zz()
                );            
            }
            
            Info << "moment of inertia (diag tensor): " << momentOfInertia_[i] << endl;
        }    
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantMoleculeProperties::~constantMoleculeProperties()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool constantMoleculeProperties::linearMoleculeTest(const label& idI) const
{
    label nSites = siteNames_[idI].size();
    
    if (nSites == 2)
    {
        return true;
    }
    
    vector refDir = siteRefPositions_[idI][1] - siteRefPositions_[idI][0];

    refDir /= mag(refDir);

    for
    (
        label i = 2;
        i < nSites;
        i++
    )
    {
        vector dir =  siteRefPositions_[idI][i] -  siteRefPositions_[idI][i-1];

        dir /= mag(dir);

        if (mag(refDir & dir) < 1 - SMALL)
        {
            return false;
        }
    }

    return true;
}

const List<diagTensor>& constantMoleculeProperties::momentOfInertia() const
{
    return momentOfInertia_;
}

const diagTensor& constantMoleculeProperties::momentOfInertia(const label& id) const
{
    return momentOfInertia_[id];
}

bool constantMoleculeProperties::linearMolecule(const label& id) const
{
    return ((momentOfInertia_[id].xx() < 0) && (momentOfInertia_[id].yy() > 0));
}

bool constantMoleculeProperties::pointMolecule(const label& id) const
{
    return (momentOfInertia_[id].zz() < 0);
}

label constantMoleculeProperties::degreesOfFreedom(const label& id) const
{
    if (linearMolecule(id))
    {
        return 5;
    }
    else if (pointMolecule(id))
    {
        return 3;
    }
    else
    {
        return 6;
    }
}

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

const List<scalar>& constantMoleculeProperties::mass() const
{
    return masses_;
}    

const scalar& constantMoleculeProperties::mass(const label& id) const
{
    return masses_[id];    
}

label constantMoleculeProperties::nSites(const word& idName) const
{
    return siteNames(idName).size();
}

label constantMoleculeProperties::nSites(const label& id) const
{
    return siteNames_[id].size();
}

const  List<List<word> >& constantMoleculeProperties::siteNames() const
{
    return siteNames_;
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

const List<List<vector> >& constantMoleculeProperties::siteRefPositions() const
{
    return siteRefPositions_;
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

        
const List<List<label> >& constantMoleculeProperties::pairPotNamesToPairPotSitesList() const
{
    return pairPotNamesToPairPotSitesList_;
}

const List<List<label> >& constantMoleculeProperties::pairPotNamesToSites() const
{
    return pairPotNamesToSites_;
}

const List<List<label> >& constantMoleculeProperties::chargePotNamesToChargePotSitesList() const
{
    return chargePotNamesToChargePotSitesList_;
}

const List<List<label> >& constantMoleculeProperties::chargePotNamesToSites() const
{
    return chargePotNamesToSites_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
