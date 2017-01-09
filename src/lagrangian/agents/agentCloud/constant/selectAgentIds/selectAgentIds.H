/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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


Description

SourceFiles
    selectAgentIdsI.H
    selectAgentIds.C
    selectAgentIdsIO.C

\*---------------------------------------------------------------------------*/

#ifndef selectAgentIds_H
#define selectAgentIds_H

#include "polyMesh.H"
#include "agentProperties.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class selectAgentIds Declaration
\*---------------------------------------------------------------------------*/

class selectAgentIds
{
    // Private data

        List<label> agentIds_;

        List<word> agentIdNames_;

    // Private Member Functions

public:

    // Static data members

    // Constructors

        //- Construct null
        selectAgentIds();

        selectAgentIds
        (
            const agentProperties& cP,
            const dictionary& dict
        );

        selectAgentIds
        (
            const agentProperties& cP,
            const dictionary& dict,
            const word& molIdsHeader
        );
        
    // Selectors

    // Destructor

        ~selectAgentIds();

    // Member Functions

        // Access

            inline List<label> agentIds();
            inline List<word> agentIdNames();            
//             inline labelList& molIds();


        // Check

        // Edit

        // Write


    // Member Operators

//         void operator=(const selectAgentIds&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const selectAgentIds&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "selectAgentIdsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
