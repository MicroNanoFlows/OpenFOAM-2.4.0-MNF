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
    pairPotentials

Description


\*---------------------------------------------------------------------------*/

#ifndef pairPotentials_H
#define pairPotentials_H

#include "vector.H"
#include "List.H"
#include "IOdictionary.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "pairPotentialModel.H"
#include "reducedUnits.H"
#include "constantMoleculeProperties.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pairPotentials Declaration
\*---------------------------------------------------------------------------*/

class pairPotentials
{

    // Private data

        const polyMesh& mesh_;

        const constantMoleculeProperties& cP_;
        
        IOdictionary potentialsDict_;

        // maximum rCut
        scalar rCut_;
//         scalar rCutSqr_;
        
        // these are just pair potentials
        
        PtrList<entry> pairPotentialsList_;
        List<word> pairNames_;
        List<label> pairIds_;
        List< autoPtr<pairPotentialModel> > pairPotentials_;

        
        // linking pairpotids to potential storage
        // outer list size of cP_.()
        // inner list size of cP_.()
        // store index of pairPotentials_
        List<List<label> > pairPotLinks_;


        
        // these are electrostatic potentials
        
        

        
        
        
        // private member functions 
        
            void testPairPotentials();
            
        
public:

    // Constructors
    
    
        //- Standard constructor
        pairPotentials
        (
            const polyMesh& mesh,
            const constantMoleculeProperties& cP,
            const reducedUnits& redUnits
        );


        ~pairPotentials();

    // Member Functions

        scalar maxRCut();
        
        inline List< autoPtr<pairPotentialModel> >& pairPots();
        inline const List< autoPtr<pairPotentialModel> >& pairPots() const;
        inline const List<word>& pairNames() const;
        inline const List<label>& pairIds() const;
        
        scalar pairForce
        (
            const label& pairSiteA,
            const label& pairSiteB,
            const scalar& r
        );
        
        scalar pairEnergy
        (
            const label& pairSiteA,
            const label& pairSiteB,
            const scalar& r
        );

        bool rCutSqr
        (
            const label& pairSiteA,
            const label& pairSiteB,
            const scalar& r
        );         
        
    // Ostream Operator

//         friend Ostream& operator<<(Ostream&, const pairPotentials&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pairPotentialsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
