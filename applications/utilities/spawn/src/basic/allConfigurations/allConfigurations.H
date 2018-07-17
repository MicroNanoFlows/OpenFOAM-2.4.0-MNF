/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description



\*---------------------------------------------------------------------------*/

#ifndef allConfigurations_H
#define allConfigurations_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "configuration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class allConfigurations Declaration
\*---------------------------------------------------------------------------*/

class allConfigurations
{
    // Private data

        //- The entire dictionary (containing multiple subDictionaries)
        const IOdictionary& dict_;
        
        PtrList<entry> configurationList_;
        List<label> ids_;
        List<autoPtr<configuration> > configurations_;

public:

    // Constructors

        //- Null Constructor
        allConfigurations
        (
//             const polyMesh& mesh,
            const IOdictionary& dict
        );

        //- Constructor for gnemdFOAM
        allConfigurations
        (
//             const polyMesh& mesh,
            molecules& cloud,
            const IOdictionary& dict
        );


    // Member Functions

        //- generates the molecules
        void spawn();

        // Access

        inline List<autoPtr<configuration> >& configurations();
        
        inline const List<autoPtr<configuration> >& configurations() const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "allConfigurationsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
