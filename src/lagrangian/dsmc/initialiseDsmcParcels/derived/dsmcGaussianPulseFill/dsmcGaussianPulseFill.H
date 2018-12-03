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
    dsmcGaussianPulseFill

Description


SourceFiles
    dsmcGaussianPulseFillI.H
    dsmcGaussianPulseFill.C
    dsmcGaussianPulseFillIO.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcGaussianPulseFill_H
#define dsmcGaussianPulseFill_H

#include "dsmcAllConfigurations.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcGaussianPulseFill Declaration
\*---------------------------------------------------------------------------*/

class dsmcGaussianPulseFill
:
	public dsmcConfiguration
{

private:

    // Private data

        //- subDictionary containing the properties
//         dictionary propsDict_;


    // Private Member Functions

    scalar gaussianAmplitude_;
    scalar sigma_;
    scalar startPoint_;


public:

    //- Runtime type information
    TypeName("dsmcGaussianPulseFill");

    // Constructors

        //- Construct from components
        dsmcGaussianPulseFill
		(
            dsmcCloud& cloud,
            const dictionary& dict
		);

    // Destructor

 		~dsmcGaussianPulseFill();


    // Member Functions

        void setInitialConfiguration();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
