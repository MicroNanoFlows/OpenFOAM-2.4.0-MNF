/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::abInitio

Description
    Variable Hard Sphere BinaryCollision Model

\*---------------------------------------------------------------------------*/

#ifndef abInitio_H
#define abInitio_H

#include "BinaryCollisionModel.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class abInitio Declaration
\*---------------------------------------------------------------------------*/


class abInitio
:
    public BinaryCollisionModel
{
    
private:
    
    // Private data

        //- The coefficients dictionary
        const dictionary coeffDict_;

        //- Scattering angle matrix
        List<scalarField> xi_;
        
        scalarField sigmaT_;
        
        scalarField g_;
        
        scalar G_;
        
        label nRows_;
        
        label nColumns_;
        
        word xiTableFileName_;

        

public:

    //- Runtime type information
    TypeName("abInitio");


    // Constructors

        //- Construct from dictionary
        abInitio
        (
            const dictionary& dict,
            dsmcCloud& cloud
        );


    //- Destructor
    ~abInitio();


    // Member Functions

        //- Flag to indicate whether model activates collision model
        bool active() const;

        //- Return the collision cross section * relative velocity product
        scalar sigmaTcR
        (
            const dsmcParcel& pP,
            const dsmcParcel& pQ
        ) const;

        //- Apply collision
        void collide
        (
            dsmcParcel& pP,
            dsmcParcel& pQ,
            label& cellI
        );

        const dictionary& coeffDict() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "abInitio.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
