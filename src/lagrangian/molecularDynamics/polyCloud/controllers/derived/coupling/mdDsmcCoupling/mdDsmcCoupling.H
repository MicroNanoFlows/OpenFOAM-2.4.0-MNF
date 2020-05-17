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
    mdDsmcCoupling

Description


SourceFiles
    mdDsmcCoupling.C

\*---------------------------------------------------------------------------*/

#ifndef mdDsmcCoupling_H
#define mdDsmcCoupling_H

#include "polyCouplingController.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mdDsmcCoupling Declaration
\*---------------------------------------------------------------------------*/

class mdDsmcCoupling
:
    public polyCouplingController
{

    struct coupledMolecule
    {
        word molType;
        point position;
        vector velocity;

        coupledMolecule() :
            molType(),
            position(vector::zero),
            velocity(vector::zero)
        {
        }
    };

    struct moleculeInsert
    {
        label addedIterations;
        polyMolecule* mol;

        moleculeInsert() :
            addedIterations(0),
            mol(NULL)
        {
        }
    };

    struct cplMoleculeInsert
    {
        label nmolsInserted;
        label nIts;

        cplMoleculeInsert() :
            nmolsInserted(0),
            nIts(0)
        {
        }
    };

private:

    // Private data

        //- subDictionarys containing the properties
        dictionary propsDict_;
        dictionary propsDictSend_;
        dictionary propsDictRecv_;

        List<label> molIds_;
        List<word> molNames_;

        bool output_;

#ifdef USE_MUI
        //- Calculated cell centres
        List<mui::point3d> cellCentres_;
#endif
        //- Sending properties
#ifdef USE_MUI
        DynamicList<mui::uniface<mui::config_3d>* > sendInterfaces_; //- 3D MUI interfaces to send through
#endif
        bool sendingRegion_;          //- Determine whether this controller will send for the fixed region
        bool sendingBound_;           //- Determine whether this controller will send for the fixed boundary
        DynamicList<Foam::word> sendInterfaceNames_;

        //- Receiving properties
#ifdef USE_MUI
        DynamicList<mui::uniface<mui::config_3d>* > recvInterfaces_; //- 3D MUI interfaces to receive from
#endif
        bool receivingRegion_;          //- Determine whether this controller will receive for the fixed region
        bool receivingBound_;           //- Determine whether this controller will receive for the fixed boundary
        DynamicList<Foam::word> recvInterfaceNames_;

        const List<word>& idList;
        const reducedUnits& rU_;
        List<List<polyMolecule*> > molHistory_;
        List<List<label> > molId_;
#ifdef USE_MUI    	
    	mui::chrono_sampler_exact3d* chrono_sampler;
#endif
    	bool couplingRegion_;
        vector couplingRegionMin_;
        vector couplingRegionMax_;
        bool couplingBounds_;
        vector couplingBoundMin_;
        vector couplingBoundMax_;
        vector couplingBoundNorm_;
        vector couplingBoundZeroThick_;
        scalar overlapEnergyLimit_;
        label overlapIterations_;
        label prevMolCount_;
        label currIteration_;
        DynamicList<label> intersectingCells_;
        DynamicList<coupledMolecule> molsToSend_;
        DynamicList<coupledMolecule> molsReceived_;
        scalar boundCorr_;
        point meshMin_;
        point meshMax_;
#ifdef USE_MUI
        List<std::vector<mui::point3d> > rcvPoints_;
#endif
        List<std::vector<std::string> > rcvMolType_;
        List<std::vector<label> > rcvMolId_;
        List<std::vector<scalar> > rcvVelX_;
        List<std::vector<scalar> > rcvVelY_;
        List<std::vector<scalar> > rcvVelZ_;
        scalar initTemperature_;
        scalar initKe_;
        scalar initTemperatureDSMC_;
        scalar initKeDSMC_;
        label nparcsRcv_;

    // Private Member Functions

    moleculeInsert insertMolecule
	(
	    point& position,
	    const label& id,
	    const bool ghost,
	    vector& velocity
	);

    polyMolecule* checkForOverlaps(polyMolecule* newMol, const scalar& potEnergyLimit);
    scalar calcTemperature();
    scalar calcAvgLinearKe();
    void scaleVelocity(scalar scaleValue);

public:

    //- Runtime type information
    TypeName("mdDsmcCoupling");

    // Constructors

        //- Construct from components
        mdDsmcCoupling
        (
            Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict,
            couplingInterface2d &twoDInterfaces,
            couplingInterface3d &threeDInterfaces
        );

    // Destructor

        ~mdDsmcCoupling();


    // Member Functions

        //- create an initial configuration
        bool initialConfiguration(label stage);

        //- run at start of time step before molecules move
        bool controlAfterMove(label stage);

        //- run after forces are calculated
        void controlAfterForces();

        //- receive molecules that fall into the coupled region
        bool receiveCoupledRegion(bool init);

        //- send accelerations for ghost molecules in coupled region
        void sendCoupledRegionAcc();

        //- Collate and delete molecules that have crossed a coupling boundary
        bool findCoupledMolecules();

        //- send molecules that have crossed a coupling boundary
        label sendCoupledMolecules();

        //- receive molecules that have crossed a coupling boundary
        label receiveCoupledParcels();

        //- Insert parcels that have been received through a coupling boundary
        cplMoleculeInsert insertCoupledMolecules();

        //- output data (optional)
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- update properties from a modified dictionary (optional)
        void updateProperties(const dictionary&);

        //- Create a coupling barrier at current coupling time for all sending interfaces
        void barrier();

        //- Create a coupling barrier at specific coupling iteration for all sending interfaces
        void barrier(label iteration);

        //- Create a coupling barrier at specific coupling iteration for specific sending interface
        void barrier(label iteration, label interface);

        //- Forget the time frame received at current coupling iteration for all interfaces
        void forget(bool forget);

        //- Forget the time frame received at specific iteration for all interfaces
        void forget(label iteration, bool forget);

        //- Forget the time frame received at specific iteration for specific interfaces
        void forget(label iteration, label interface, bool forget);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
