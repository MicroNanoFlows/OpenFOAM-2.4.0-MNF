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
    dsmcCouplingController

Description
    Basic/abstract class of a coupling controller


SourceFiles
    dsmcCouplingController.C
    dsmcCouplingControllerIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcCouplingController_H
#define dsmcCouplingController_H

#include "coupling2d.H"
#include "coupling3d.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "Random.H"
#include "dsmcParcel.H"
#include "timeData.H"
#include "selectIds.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcCouplingController Declaration
\*---------------------------------------------------------------------------*/

class dsmcCouplingController
{

protected:

    // Protected data


        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        Random& rndGen_;

        timeData time_;

        // on/off
        bool control_;

        //- writing out properties
        bool writeInTimeDir_;
        bool writeInCase_;

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("dsmcCouplingController");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            dsmcCouplingController,
            dictionary,
            (
                Time& t,
                dsmcCloud& cloud,
                const dictionary& dict,
                couplingInterface2d& twoDInterfaces,
                couplingInterface3d& threeDInterfaces
            ),
            (t, cloud, dict, twoDInterfaces, threeDInterfaces)
        );
   
    // Constructors

        //- Construct from components
        dsmcCouplingController
        (
			Time& t,
            dsmcCloud& cloud,
            const dictionary& dict,
            couplingInterface2d& twoDInterfaces,
            couplingInterface3d& threeDInterfaces
        );


    // Selectors

        static autoPtr<dsmcCouplingController> New
        (
			Time& t,
            dsmcCloud& cloud,
            const dictionary& dict,
            couplingInterface2d& twoDInterfaces,
            couplingInterface3d& threeDInterfaces
        );


    // Destructor

        virtual ~dsmcCouplingController();


    // Member Functions

        //- create an initial configuration
        virtual void initialConfiguration(label stage) = 0;

        //- run after move and before collisions
        virtual void controlParcelsBeforeCollisions(label stage) = 0;

        //- run after move and collisions
        virtual void controlParcelsAfterCollisions(label stage) = 0;

        //- send parcels that fall into the MD coupled region to coupled mdFoamPlus
        virtual void sendCoupledRegion(bool init) = 0;

        //- send parcels to mdFoamPlus that have crossed a coupling boundary
        virtual void sendCoupledParcels() = 0;

        //- receive molecules from mdFoamPlus that have crossed a coupling boundary
        virtual bool receiveCoupledMolecules() = 0;

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePaths
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
        virtual void forget(bool forget) = 0;

        //- Forget the time frame received at specific time for all interfaces
        virtual void forget(label iteration, bool forget) = 0;

        //- Forget the time frame received at specific time for specific interfaces
        virtual void forget(label iteration, label interface, bool forget) = 0;

        //- Reset ghosted status of cloud
        virtual void resetGhostedStatus() = 0;

    // Access

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
