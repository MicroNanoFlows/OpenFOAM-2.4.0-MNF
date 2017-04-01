/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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
    Foam::polyTracker

Description



SourceFiles
    polyTrackerI.H
    polyTracker.C
    polyTrackerIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyTracker_H
#define polyTracker_H

#include "scalar.H"
#include "IOdictionary.H"
#include "polyMesh.H"
#include "Time.H"
#include "polyMolecule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class polyTracker Declaration
\*---------------------------------------------------------------------------*/

class polyTracker
{
    // Private data

        polyMoleculeCloud& cloud_;
    
        DynamicList<polyMolecule*> molsOnProc_;
        
        DynamicList<label> molTrackingNumbers_;
        

        // tool variables
        
//         label trackingNumber_;
        
        label id_;
        

    // Private Member Functions



public:
    
        friend class polyMolecule;


    // Constructors


        //- Construct from components
        polyTracker
        (
            polyMoleculeCloud& cloud
        );


    // Destructor

        ~polyTracker();


    // Member Functions

            void track();
        
        // Access
            
            // this is a two step function
            
            inline label size() const;
            
            
//             inline const DynamicList<label>& indices() const;
            
            
            // call this first (to avoid errors)
            inline bool isTrackingNumberAvailable
            (
                const label& trackingNumber
            );

            // then call this function to get the polyMolecule
            inline polyMolecule* getMol();

            inline const DynamicList<label>& molTrackingNumbers() const;            
            
            inline DynamicList<polyMolecule*>& molsOnProc();
            inline const DynamicList<polyMolecule*>& molsOnProc() const;

            inline label noOfMolsOnProc(); 
            
    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const polyTracker&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyTrackerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
