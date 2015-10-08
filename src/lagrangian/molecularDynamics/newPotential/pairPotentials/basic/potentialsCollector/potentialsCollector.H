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
    potentialsCollector

Description


\*---------------------------------------------------------------------------*/

#ifndef potentialsCollector_H
#define potentialsCollector_H

#include "vector.H"
#include "List.H"
#include "IOdictionary.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "pairPotentialModel.H"
#include "reducedUnits.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class potentialsCollector Declaration
\*---------------------------------------------------------------------------*/

class potentialsCollector
{

    // Private data

        const polyMesh& mesh_;

        IOdictionary potentialsDict_;

        // these are just pair potentials
        
        PtrList<entry> pairPotentialsList_;
        List<word> pairNames_;
        List<label> pairIds_;
        List< autoPtr<pairPotentialModel> > pairPotentials_;

        
        // any different catagories of potentials can go here
        
        
        
        // private member functions 
        
//             void testPairPotentials();
        
            
        
public:

    // Constructors

        //- Null Constructor (mdInitialise)
        potentialsCollector
        (
            const polyMesh& mesh,
            const reducedUnits& redUnits
        );


        ~potentialsCollector();

    // Member Functions

        inline List< autoPtr<pairPotentialModel> >& pairPotentials();
        inline const List< autoPtr<pairPotentialModel> >& pairPotentials() const;
        inline const List<word>& pairNames() const;
        inline const List<label>& pairIds() const;
        
        
    // Ostream Operator

//         friend Ostream& operator<<(Ostream&, const potentialsCollector&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "potentialsCollectorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
