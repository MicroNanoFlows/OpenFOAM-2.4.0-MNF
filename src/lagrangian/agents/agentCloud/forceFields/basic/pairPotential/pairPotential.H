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
    pairPotential

Description
    Basic/abstract class of a state controller


SourceFiles
    pairPotentialI.H
    pairPotential.C
    pairPotentialIO.C

\*---------------------------------------------------------------------------*/

#ifndef pairPotential_H
#define pairPotential_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "agent.H"
#include "writeTimeData.H"
#include "selectAgentIds.H"
#include "reducedUnits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pairPotential Declaration
\*---------------------------------------------------------------------------*/

class pairPotential
{

protected:

    // Protected data

        const fvMesh& mesh_;

        agentCloud& cloud_;

        const reducedUnits& rU_;        
        
        // ids header of dictionary e.g. A-B 
        word name_;        
        
        // split ids [0] = A; [1] = B
        List<word> idList_;
        
        scalar rCut_;
        scalar rCutSqr_;
        scalar rMin_;
        scalar dr_;           
        
//         bool writeTable_;
        
        //- subDictionary containing the properties
//         dictionary controllerDict_;

//         Time& time_; 



    // Private Member Functions

public:

    //- Runtime type information
    TypeName("pairPotential");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            pairPotential,
            dictionary,
            (
                agentCloud& cloud,
                const word& name,
                const dictionary& dict
            ),
            (cloud, name, dict)
        );
   
    // Constructors

        //- Construct from components
        pairPotential
        (
            agentCloud& cloud,
            const word& name,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<pairPotential> New
        (
            agentCloud& cloud,
            const word& name,
            const dictionary& dict
        );


    // Destructor

        virtual ~pairPotential();


    // Member Functions


        //- create an initial configuration
        virtual void initialConfiguration() = 0;
        
        virtual scalar energy(const scalar& r) = 0;
        
        virtual scalar force(const scalar& r) = 0;
        
        virtual void pairPotentialFunction
        (
            agent* molI,
            agent* molJ,
            const scalar& r,
            scalar& energy,
            scalar& force
        ) = 0;        

        inline const scalar& rCutSqr() const;
        
        const List<word>& idList() const;
        
        void writeTables(const fileName& pathName);
        
        

/*
        //- control molecules at various stages of time-step
        virtual void controlBeforeVelocityI() = 0;
        
        virtual void controlBeforeMove() = 0;
        
        virtual void controlBeforeForces() = 0;
        
        virtual void controlDuringForces
        (
            agent* molI,
            agent* molJ
        ) = 0;
        
        virtual void controlAfterForces() = 0;
        
        virtual void controlAfterVelocityII() = 0;

        //- calculate any required properties
        virtual void calculateProperties() = 0;
        
        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;*/

        //- update properties from a modified dictionary
//         virtual void updateProperties(const dictionary&) = 0;

//         void updateStateControllerProperties(const dictionary&);

        // Access

//             //- return the control zone cells
//             const labelList& controlZone() const;
// 
//             //- return the control zone name
//             const word& regionName() const;
// 
//             const bool& controlInterForces() const;
//             bool& controlInterForces();
// 
//             const bool& writeInTimeDir() const;
//             const bool& writeInCase() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pairPotentialI.H"

#endif

// ************************************************************************* //
