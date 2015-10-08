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
    REQUIRES A RENAME

    Collect information on the molecules being simulated
    Read in from constant/moleculeProperties
    
    tetherIdList can be removed from this file 

    use SI units

SourceFiles
    constantMoleculeProperties.C

\*---------------------------------------------------------------------------*/

#ifndef constantMoleculeProperties_H
#define constantMoleculeProperties_H

#include "Time.H"
#include "polyMesh.H"

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
        
        // number of molecules types 
        label N_;
        
        // ordered list of molecule names
        List<word> names_;
        
        // outer list = molecule types, linked to cloud (e.g. polyMoleculeCloud)
        List<word> cloudTypes_;

        // outer list = molecule types, inner list sites corresponding to  
        List<List<word> > siteNames_;

        List<List<word> > pairPotNames_;        

        List<List<vector> > siteRefPositions_;

        List<List<scalar> > siteMasses_;

        List<List<scalar> > siteCharges_;
        
        label nPairPotIds_;

        // ordered lists
        // full list of sites (unique list)
        List<word> siteIdList_;
        
        // just those that are pair potential (unique list)
        List<word> potSiteIdList_;
        
        
public:

    // Constructors

        //- constructor
        constantMoleculeProperties(const polyMesh& mesh);

    // Destructor

        ~constantMoleculeProperties();


    // Member functions

        // no of id types
        const label& nTypes() const;
        
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
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
