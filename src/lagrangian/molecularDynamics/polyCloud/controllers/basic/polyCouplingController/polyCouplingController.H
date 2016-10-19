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
    polyStateController

Description
    Basic/abstract class of a state controller


SourceFiles
    polyCouplingController.C

\*---------------------------------------------------------------------------*/

#ifndef polyCouplingController_H
#define polyCouplingController_H

#include "coupling1d.H"
#include "coupling2d.H"
#include "coupling3d.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "polyMolecule.H"
#include "writeTimeData.H"
#include "selectIds.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyCouplingController Declaration
\*---------------------------------------------------------------------------*/

class polyCouplingController
{

protected:

    // Protected data

        const fvMesh& mesh_;

        polyMoleculeCloud& molCloud_;

        Time& time_; 

        //- name of control zone
        word regionName_;
        label regionId_;

        // on/off
        bool control_;

        // expensive control within intermolecular force computation
        bool controlInterForces_;

        bool writeInTimeDir_;
        bool writeInCase_;

        //- MUI coupling interfaces
		couplingInterface1d oneDInterfaces_;
		couplingInterface2d twoDInterfaces_;
		couplingInterface3d threeDInterfaces_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("polyCouplingController");

    // Declare runtime constructor selection table
    	declareRunTimeSelectionTable
        (
            autoPtr,
			polyCouplingController,
            dictionary,
            (
                Time& t,
                polyMoleculeCloud& molCloud,
                const dictionary& dict,
				couplingInterface1d &oneDInterfaces,
				couplingInterface2d &twoDInterfaces,
				couplingInterface3d &threeDInterfaces
            ),
            (t, molCloud, dict, oneDInterfaces, twoDInterfaces, threeDInterfaces)
        );
   
    // Constructors

        //- Construct from components
        polyCouplingController
        (
            Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict,
			couplingInterface1d &oneDInterfaces,
			couplingInterface2d &twoDInterfaces,
			couplingInterface3d &threeDInterfaces
        );


    // Selectors

        static autoPtr<polyCouplingController> New
        (
            Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict,
			couplingInterface1d &oneDInterfaces,
			couplingInterface2d &twoDInterfaces,
			couplingInterface3d &threeDInterfaces
        );


    // Destructor

        virtual ~polyCouplingController();


    // Member Functions


        //- create an initial configuration
        virtual void initialConfiguration() = 0;

        //- receive any coupled data to the MUI interface (start of time-step)
		virtual void receiveCoupling() = 0;

		//- send any coupled data to the MUI interface (end of time-step)
		virtual void sendCoupling() = 0;
        
        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        //- update properties from a modified dictionary
        virtual void updateProperties(const dictionary&) = 0;

        void updateCouplingControllerProperties(const dictionary&);

	// Access

		//- return the control zone cells
		const labelList& controlZone() const;

		//- return the control zone name
		const word& regionName() const;

		const bool& controlInterForces() const;
		bool& controlInterForces();

		const bool& writeInTimeDir() const;
		const bool& writeInCase() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
