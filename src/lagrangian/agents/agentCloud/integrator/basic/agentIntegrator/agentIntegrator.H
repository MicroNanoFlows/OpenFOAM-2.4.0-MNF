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
    agentIntegrator

Description
    Basic/abstract class of a state controller


SourceFiles
    agentIntegratorI.H
    agentIntegrator.C
    agentIntegratorIO.C

\*---------------------------------------------------------------------------*/

#ifndef agentIntegrator_H
#define agentIntegrator_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "agent.H"
#include "writeTimeData.H"
#include "selectIds.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agentIntegrator Declaration
\*---------------------------------------------------------------------------*/

class agentIntegrator
{

protected:

    // Protected data

        const fvMesh& mesh_;

        agentCloud& cloud_;

        Time& time_; 

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("agentIntegrator");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            agentIntegrator,
            dictionary,
            (
                Time& t,
                agentCloud& cloud,
                const dictionary& dict
            ),
            (t, cloud, dict)
        );
   
    // Constructors

        //- Construct from components
        agentIntegrator
        (
            Time& t,
            agentCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<agentIntegrator> New
        (
            Time& t,
            agentCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~agentIntegrator();


    // Member Functions


        //- create an initial configuration
        virtual void init() = 0;                
        
        //- evolve MD simulation
        virtual void evolve() = 0;
        
        virtual void clearLagrangianFields() = 0; 

        virtual void calculateForce() = 0; 
        
        virtual void updateAcceleration() = 0;         
    
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
