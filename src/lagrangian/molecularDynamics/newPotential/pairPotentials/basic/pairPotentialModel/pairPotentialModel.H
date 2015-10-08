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
    pairPotentialModel

Description

SourceFiles
    pairPotentialModelI.H
    pairPotentialModel.C
    pairPotentialModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef pairPotentialModel_H
#define pairPotentialModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "writeTimeData.H"
#include "polyMesh.H"
#include "mathematicalConstants.H"
#include "reducedUnits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Forward declaration of classes
class energyScalingFunction;

/*---------------------------------------------------------------------------*\
                           Class pairPotentialModel Declaration
\*---------------------------------------------------------------------------*/

class pairPotentialModel
{

protected:

    // Protected data

        const polyMesh& mesh_;

        const reducedUnits& rU_;
        
        scalar rCut_;
        scalar rCutSqr_;
        
        scalar rMin_;
        
        scalar dr_;        
        
//         word energyScalingFunction_;


        scalar fMin_;
        scalar energyMin_;        
       
        List<scalar> forceLookup_;
        List<scalar> energyLookup_;
        
        mutable energyScalingFunction* esfPtr_;

        bool writeTables_;        

        
        
    // Private Member Functions
        
        void scaleEnergy
        (
            scalar& e,
            const scalar r,
            const reducedUnits& rU
        ) const;        

public:

    //- Runtime type information
    TypeName("pairPotentialModel");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            pairPotentialModel,
            dictionary,
            (
                const polyMesh& mesh,
                const reducedUnits& redUnits,
                const dictionary& dict
            ),
            (mesh, redUnits, dict)
        );
   
    // Constructors

        //- Construct from components
        pairPotentialModel
        (
            const polyMesh& mesh,
            const reducedUnits& redUnits,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<pairPotentialModel> New
        (
            const polyMesh& mesh,
            const reducedUnits& redUnits,
            const dictionary& dict
        );



    // Destructor

        virtual ~pairPotentialModel();


    // Member Functions
        
        void setLookupTables();

        inline scalar rMin() const;

        inline scalar dr() const;

        inline scalar rCut() const;

        inline scalar rCutSqr() const;

        scalar energy (const scalar r) const;

        scalar force (const scalar r) const;

        List<Pair<scalar> > energyTable() const;

        List<Pair<scalar> > forceTable() const;

        inline bool writeTables() const;

        virtual scalar unscaledEnergy(const scalar r) const = 0;

        virtual const dictionary& dict() const = 0;

        scalar scaledEnergy
        (
            const scalar r
        ) const;

        scalar energyDerivative
        (
            const scalar r,
            const bool scaledEnergyDerivative = true
        ) const;

//         const dictionary& pairPotentialProperties() const
//         {
//             return pairPotentialProperties_;
//         }

        bool writeEnergyAndForceTables(Ostream& os) const;

        //- Read pairPotential dictionary
//         virtual bool read
//         (
//             const dictionary& pairPotentialProperties
//         ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "pairPotentialModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
