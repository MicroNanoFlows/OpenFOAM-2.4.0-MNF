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
    helbingVelocityDependentBAnisotropy

Description
    Model that applies a gravitational force on all molecules.

SourceFiles
    helbingVelocityDependentBAnisotropyI.H
    helbingVelocityDependentBAnisotropy.C
    helbingVelocityDependentBAnisotropyIO.C

\*---------------------------------------------------------------------------*/

#ifndef helbingVelocityDependentBAnisotropy_H
#define helbingVelocityDependentBAnisotropy_H

#include "pairPotential.H"
#include "vector.H"
#include "polyMesh.H"
#include "agentCloud.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class helbingVelocityDependentBAnisotropy Declaration
\*---------------------------------------------------------------------------*/

class helbingVelocityDependentBAnisotropy
:
	public pairPotential
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        scalar A_;
        scalar Bmax_;
	scalar Bmin_;
        scalar k_;
        scalar kappa_;
	scalar lambda_;


        
    // Private Member Functions

 

public:

    //- Runtime type information
    TypeName("helbingVelocityDependentBAnisotropy");

    // Constructors

        //- Construct from components
        helbingVelocityDependentBAnisotropy
        (
            agentCloud& cloud,
            const word& name,
            const dictionary& dict
        );

    // Destructor

        ~helbingVelocityDependentBAnisotropy();


    // Member Functions
    // Member Functions
        
        //- create an initial configuration
        void initialConfiguration();
        
        scalar energy(const scalar& r);
        
        scalar force(const scalar& r);
        
        void pairPotentialFunction
        (
            agent* molI,
            agent* molJ,
            const scalar& r,
            scalar& energy,
            vector& force
        );        
        
//         void writeTables(const fileName& pathName);

        // Access
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
