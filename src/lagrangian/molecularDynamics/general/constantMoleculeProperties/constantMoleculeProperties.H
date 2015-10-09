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

Class
    constantMoleculeProperties

Description


SourceFiles
    constantMoleculeProperties.C

\*---------------------------------------------------------------------------*/

#ifndef constantMoleculeProperties_H
#define constantMoleculeProperties_H

#include "Time.H"
#include "polyMesh.H"
#include "reducedUnits.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class constantMoleculeProperties Declaration
\*---------------------------------------------------------------------------*/

class constantMoleculeProperties
{
    // Private data

        const polyMesh& mesh_;
        
        const reducedUnits& rU_;
        
        // number of molecules types 
        label N_;
        
        // ordered list of molecule names
        List<word> names_;
        
        // outer list = molecule types, linked to cloud (e.g. polyMoleculeCloud)
        List<word> cloudTypes_;

        // outer list = molecule types
        // inner list = all sites in molecule 
        // entry = name
        List<List<word> > siteNames_;

        // outer list = molecule types
        // inner list = just pair-pot sites in molecule 
        // entry = name
        List<List<word> > pairPotNames_;        

        // outer list = molecule types
        // inner list = just charge sites in molecule // duplicates allowed
        // entry = name
        List<List<word> > chargeNames_;         
        
        // outer list = molecule types
        // inner list = all sites in molecule 
        // entry = position of site
        List<List<vector> > siteRefPositions_;

        // outer list = molecule types
        // inner list = all sites in molecule 
        // entry = site mass
        List<List<scalar> > siteMasses_;
        
        // outer list = molecule types
        // inner list = all sites in molecule 
        // entry = site charge
        List<List<scalar> > siteCharges_;
        
        label nSites_;
        label nPairPotSites_;
        label nChargeSites_;        
        
        // ordered unqiue lists (A,B,C, ..., N)
       
        // full list of sites (unique list)
        // Outer list = total number of sites in system
        // Entry = name
        List<word> siteIdList_;
        
        // just those that are pair potential (unique list)
        // Outer list = total number of pair pot sites in system
        // Entry = name
        List<word> pairPotSiteIdList_;
        
        // just those sites that are charged (unique list)
        // Outer list = total number of non-zero charged sites in system
        // Entry = name        // duplicates not allowed
        List<word> chargeSiteIdList_;


        // linked lists
        
        // link between pairPotNames_ and pairPotSiteIdList_
        // outer list = molecule types
        // inner list = linked to pairPotNames_
        // each entry is an index to pairPotSiteIdList_
        List< List<label> > pairPotNamesToFullSites_;

        // link between pairPotNames_ and siteNames_
        // outer list = molecule types
        // inner list = linked to pairPotNames_
        // each entry is an index to siteNames_
        List< List<label> > pairPotNamesToSites_;   
        
        
        // link between chargeNames_ and chargeSiteIdList_
        // outer list = molecule types
        // inner list = linked to chargeNames_
        // each entry is an index to chargeSiteIdList_
        List< List<label> > chargePotNamesToFullSites_;        
        
        
public:

    // Constructors

        //- constructor
        constantMoleculeProperties
        (
            const polyMesh& mesh,
            const reducedUnits& rU
        );

    // Destructor

        ~constantMoleculeProperties();


    // Member functions

        // no of id types
        const label& nMolTypes() const;
        
        const List<word>& molIds() const;
        const List<word>& cloudTypes() const;
        const List<word>& siteNames(const word& idName) const;
        const List<word>& siteNames(const label& id) const;

        const List<word>& pairPotNames(const word& idName) const;
        const List<word>& pairPotNames(const label& id) const;

        const List<vector>& siteRefPositions(const word& idName) const;
        const List<vector>& siteRefPositions(const label& id) const;

        const List<scalar>& siteMasses(const word& idName) const;
        const List<scalar>& siteMasses(const label& id) const;
        
        const List<scalar>& siteCharges(const word& idName) const;
        const List<scalar>& siteCharges(const label& id) const;  
        
        
        const label& nSiteTypes() const;
        const label& nPairPotTypes() const;
        const label& nChargeTypes() const;
        
        // unique lists of sites
        const List<word>& siteIds() const;
        const List<word>& pairPotSiteIdList() const;
        const List<word>& chargeSiteIdList() const;
        
        
        const List<List<label> >& pairPotNamesToFullSites() const;
        const List<List<label> >& pairPotNamesToSites() const;        
        const List<List<label> >& chargePotNamesToFullSites() const;        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
