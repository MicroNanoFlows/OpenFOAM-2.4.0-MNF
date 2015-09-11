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
    polyTemperatureBerendsenBins

Description
	The Berendsen thermostat rescales the kinetic temperature over variable time periods rather 
    than instantaneously.

SourceFiles
    polyTemperatureBerendsenBinsI.H
    polyTemperatureBerendsenBins.C
    polyTemperatureBerendsenBinsIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyTemperatureBerendsenBins_H
#define polyTemperatureBerendsenBins_H

#include "polyStateController.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "binModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyTemperatureBerendsenBins Declaration
\*---------------------------------------------------------------------------*/

class polyTemperatureBerendsenBins
:
	public polyStateController
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        //- time constant
        scalar tauT_;

        // component control
        bool componentControl_;
        bool X_;
        bool Y_;
        bool Z_;

        bool peculiar_;

        autoPtr<binModel> binModel_;

        List<label> molIds_;

        bool measureFullTemperature_;

        //- measured temperature

        scalar averagingCounter_;

        scalarField mass_;
        vectorField momentum_;
        vectorField velocity_;

        scalarField kE_;
        scalarField dof_;
        scalarField angularKe_;
        scalarField measuredT_;

		//- velocity rescaling term
		scalarField chi_;


        bool resetFieldsAtOutput_;
        bool output_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("polyTemperatureBerendsenBins");

    // Constructors

        //- Construct from components
        polyTemperatureBerendsenBins
		(
            Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyTemperatureBerendsenBins();


    // Member Functions

        void initialConfiguration();

        //- B. calculate any required properties
        void calculateProperties();

        //- C. control the mols
        void controlMols();

        //- D. output data
        void output
        (
            const fileName& fixedPathName, 
            const fileName& timePath
        );

        //- E.  update properties from a modified dictionary
        void updateProperties(const dictionary&);

        void controlMolsBeg();
        void controlMolsEnd();
        void controlBeforeForces();
        void controlDuringForces
        (
            polyMolecule* molI,
            polyMolecule* molJ
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
