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
    velocityVerlet

Description


SourceFiles
    velocityVerletI.H
    velocityVerlet.C
    velocityVerletIO.C

\*---------------------------------------------------------------------------*/

#ifndef velocityVerlet_H
#define velocityVerlet_H

#include "agentIntegrator.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "agentCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class velocityVerlet Declaration
\*---------------------------------------------------------------------------*/

class velocityVerlet
:
    public agentIntegrator
{

private:

    
        // Private data
       
    
    
public:

    //- Runtime type information
    TypeName("velocityVerlet");

    // Constructors

        //- Construct from components
        velocityVerlet
        (
            Time& t,
            agentCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~velocityVerlet();


    // Member Functions

        void init();
        
        void evolve();
        
        void updateHalfVelocity();
        
        void clearLagrangianFields(); 

        void calculateForce(); 
        
        
        // 
        
        void checkMaxVelocity();

//         void applyWillForce();        
        

        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
