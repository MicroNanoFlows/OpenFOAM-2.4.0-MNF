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
    polyCouplingPatch

Description

SourceFiles
    polyCouplingPatch.C


\*---------------------------------------------------------------------------*/

#ifndef polyCouplingPatch_H
#define polyCouplingPatch_H

#include "polyPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyCouplingPatch Declaration
\*---------------------------------------------------------------------------*/

class polyCouplingPatch
:
	public polyPatchBoundary
{

private:

    // Private data
        dictionary propsDict_;

        scalar elapsedTime_;

        scalar writeInterval_;

        scalar startTime_;

        scalar molFlux_;
        scalar massFlux_;

        scalar cumulMolFlux_;
        scalar cumulMassFlux_;

public:

    //- Runtime type information
    TypeName("polyCouplingPatch");

    // Constructors

        //- Construct from components
        polyCouplingPatch
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyCouplingPatch();


    // Member Functions

        void initialConfiguration();

        //- calculate any required properties
        void calculateProperties();

        //- control the mols
        void controlMol
        (
            polyMolecule& mol,
            polyMolecule::trackingData& td
        );

        //- output data
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
