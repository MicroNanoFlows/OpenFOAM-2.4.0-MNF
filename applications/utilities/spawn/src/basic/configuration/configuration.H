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
    configuration

Description


SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef configuration_H
#define configuration_H

#include "IOdictionary.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "molecules.H"
#include "boundsBox.H"
#include "IFstream.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class configuration Declaration
\*---------------------------------------------------------------------------*/

class configuration
{
protected:

    // Protected data

        molecules& cloud_;

        const dictionary& dict_;

        
    // Private Member Functions

        void insertAtom
        (
            const vector& position,
            const label type,
            const scalar& charge
        );
        
        void deleteAtom
        (
            label id
        );

public:

    //- Runtime type information
    TypeName("configuration");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            configuration,
            dictionary,
            (
                molecules& cloud,
                const dictionary& dict
            ),
            (cloud, dict)
        );

    // Constructors

        //- Construct from components
        configuration
        (
            molecules& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<configuration> New
        (
            molecules& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~configuration();


    // Member Functions

        virtual void spawn() = 0;

        
//         label mols() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
