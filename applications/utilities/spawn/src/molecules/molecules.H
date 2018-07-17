/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    molecules

Description
    


SourceFiles
    molecules.C

\*---------------------------------------------------------------------------*/

#ifndef molecules_H
#define molecules_H

#include "vector.H"
#include "polyMesh.H"
#include "IOdictionary.H"
#include "boundsBox.H" 
#include "IFstream.H"
#include "OFstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class molecules Declaration
\*---------------------------------------------------------------------------*/


class molecules
{
    // Private data
    
        const IOdictionary& dict_;    

        DynamicList<vector> positions_;

        DynamicList<label> types_;
        
        DynamicList<scalar> charges_;
        
        List<word> typeNames_;
        List<label> typeIds_;
        
        boundsBox bMesh_;
        

        
    // Private Member Functions
        
        void setBoundsBox
        (
            const dictionary& propsDict,
            boundsBox& bb,
            const word& name 
        );

        
        
public:

    // Constructors

        //- Construct
        molecules
        (
//             const polyMesh& mesh,
            const IOdictionary& dict         
        );


    // Destructor

        ~molecules();

    // Member Functions


        // Access

            DynamicList<vector>& positions();
            const DynamicList<vector>& positions() const;

            DynamicList<label>& types();
            const DynamicList<label>& types() const;

            DynamicList<scalar>& charges();
            const DynamicList<scalar>& charges() const;
            
            const List<word>& typeNames() const;
            const List<label>& typeIds() const;
            
            const boundsBox& bMesh() const;
            boundsBox& bMesh();
            
            
            label getType(const word& typeName);
            
            
        // Check

        // Edit

        // Write
            
            void write();
            
            label size();

    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
