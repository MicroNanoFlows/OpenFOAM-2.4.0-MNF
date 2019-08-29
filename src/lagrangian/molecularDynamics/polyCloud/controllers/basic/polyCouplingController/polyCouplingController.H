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
    polyCouplingController

Description
    Basic/abstract class of a state controller


SourceFiles
    polyCouplingController.C

\*---------------------------------------------------------------------------*/

#ifndef polyCouplingController_H
#define polyCouplingController_H

#include "coupling2d.H"
#include "coupling3d.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "polyMolecule.H"
#include "writeTimeData.H"

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
        List<word> regionNames_;
        List<label> regionIds_;

        // on/off
        bool control_;

        // expensive control within intermolecular force computation
        bool controlInterForces_;

        bool writeInTimeDir_;
        bool writeInCase_;

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
                couplingInterface2d &twoDInterfaces,
                couplingInterface3d &threeDInterfaces
            ),
            (t, molCloud, dict, twoDInterfaces, threeDInterfaces)
        );

    // Constructors

        //- Construct from components
        polyCouplingController
        (
            Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict,
            couplingInterface2d &twoDInterfaces,
            couplingInterface3d &threeDInterfaces
        );


    // Selectors

        static autoPtr<polyCouplingController> New
        (
            Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict,
            couplingInterface2d &twoDInterfaces,
            couplingInterface3d &threeDInterfaces
        );


    // Destructor

        virtual ~polyCouplingController();


    // Member Functions


        //- create an initial configuration
        virtual void initialConfiguration() = 0;

        //- run at start of MD time step before molecules move
        virtual void controlAfterMove(label stage) = 0;

        //- run after forces are calculated
        virtual void controlAfterForces() = 0;

        //- receive positions and velocities of molecules that fall into the coupled region
        virtual void receiveCoupledRegion(bool init) = 0;

        //- send forces acting on molecules in coupled region
        virtual void sendCoupledRegionForces() = 0;

        //- send molecules to dsmcFoamPlus that have crossed a coupling boundary
        virtual void sendCoupledMolecules() = 0;

        //- receive parcels from dsmcFoamPlus that have crossed a coupling boundary
        virtual void receiveCoupledParcels() = 0;

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        //- update properties from a modified dictionary
        virtual void updateProperties(const dictionary&) = 0;

        void updateCouplingControllerProperties(const dictionary&);

        //- create a coupling barrier at current coupling time for all sending interfaces
		virtual void barrier() = 0;

		//- create a coupling barrier for a specific sending interface
		virtual void barrier(label time, label interface) = 0;

		//- create a coupling barrier for all sending interfaces
		virtual void barrier(label iteration) = 0;

		//- Forget the time frame received at current coupling time for all interfaces
		virtual void forget() = 0;

		//- Forget the time frame received at specific time for all interfaces
		virtual void forget(label iteration) = 0;

		//- Forget the time frame received at specific time for specific interfaces
		virtual void forget(label iteration, label interface) = 0;

    // Access

        //- return the control zone cells
        const labelList& controlZone(label regionID) const;

        //- return the control zone names
        const List<word>& regionNames() const;

        //- return the control zone IDs
        const List<label>& regionIds() const;

        const bool& controlInterForces() const;
        bool& controlInterForces();

        const bool& writeInTimeDir() const;
        const bool& writeInCase() const;

        void rebuildCellOccupancy();
        void prepareInteractions();
        void resetTracking();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
