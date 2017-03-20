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
    helbingExponential

Description
    Model that applies a gravitational force on all molecules.

SourceFiles
    helbingExponentialI.H
    helbingExponential.C
    helbingExponentialIO.C

\*---------------------------------------------------------------------------*/

#ifndef helbingExponential_H
#define helbingExponential_H

#include "pairPotential.H"
#include "vector.H"
#include "polyMesh.H"
#include "agentCloud.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class helbingExponential Declaration
\*---------------------------------------------------------------------------*/

class helbingExponential
:
	public pairPotential
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        scalar A_;
        scalar B_;
        scalar k_;
        scalar kappa_;
        
        label option_;
                
        bool injury_;
        scalar maxForce_; // per unit length
        label agentId_;
        
        scalar timeDelay_;
        
        // anisotropic 
        scalar thetaRef1_;
        
//         scalar C_;
//         scalar thetaRef2_;
        
    // Private Member Functions

 

public:

    //- Runtime type information
    TypeName("helbingExponential");

    // Constructors

        //- Construct from components
        helbingExponential
        (
            agentCloud& cloud,
            const word& name,
            const dictionary& dict
        );

    // Destructor

        ~helbingExponential();


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
