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

#include "coupling1d.H"
#include "coupling2d.H"
#include "coupling3d.H"
#include "IOdictionary.H"
// #include "Time.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "Random.H"
#include "dsmcParcel.H"
#include "timeData.H"


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

        //- subDictionary containing the properties
        dictionary controllerDict_;

        dictionary timeDict_;

        timeData time_;

        scalar timePeriod_;
        scalar initialTime_;

        //- name of control zone
        word regionName_;
        label regionId_;

        bool control_;
        bool readStateFromFile_;

        //- set all the properties below from model if required

        // target single valued controller
        bool singleValueController_;

        scalar density_;
        vector velocity_;
        scalar temperature_;
        scalar pressure_;

        tensor strainRate_;
        vector tempGradient_;

        //- targeted field controller
        bool fieldController_;

        scalarField densities_;
        vectorField velocities_;
        scalarField temperatures_;
        scalarField pressures_;

        //- writing out properties
        bool writeInTimeDir_;
        bool writeInCase_;

        //- MUI coupling interfaces
        couplingInterface1d oneDInterfaces_;
        couplingInterface2d twoDInterfaces_;
        couplingInterface3d threeDInterfaces_;

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
                couplingInterface1d &oneDInterfaces,
                couplingInterface2d &twoDInterfaces,
                couplingInterface3d &threeDInterfaces
            ),
            (t, cloud, dict, oneDInterfaces, twoDInterfaces, threeDInterfaces)
        );
   
    // Constructors

        //- Construct from components
        dsmcCouplingController
        (
			Time& t,
            dsmcCloud& cloud,
            const dictionary& dict,
            couplingInterface1d &oneDInterfaces,
            couplingInterface2d &twoDInterfaces,
            couplingInterface3d &threeDInterfaces
        );


    // Selectors

        static autoPtr<dsmcCouplingController> New
        (
			Time& t,
            dsmcCloud& cloud,
            const dictionary& dict,
            couplingInterface1d &oneDInterfaces,
            couplingInterface2d &twoDInterfaces,
            couplingInterface3d &threeDInterfaces
        );


    // Destructor

        virtual ~dsmcCouplingController();


    // Member Functions

        //- create an initial configuration
        virtual void initialConfiguration() = 0;

        //- run at start of DSMC time step before particles move
        virtual void controlParcelsBeforeMove() = 0;

        //- run at end of DSMC time step
        virtual void calculateProperties() = 0;

        //- send parcels that fall into the MD coupled region to coupled mdFoamPlus
        virtual void sendCoupledRegion() = 0;

        //- send parcels to mdFoamPlus that have crossed a coupling boundary
        virtual void sendCoupledParcels() = 0;

        //- receive molecules from mdFoamPlus that have crossed a coupling boundary
        virtual void receiveCoupledMolecules() = 0;

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
