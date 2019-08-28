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
#include <sstream>

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
        bool sending_;          //- Determine whether this interface will send
        DynamicList<Foam::word> sendInterfaceNames_;

        //- Receiving properties
#ifdef USE_MUI
        DynamicList<mui::uniface<mui::config_3d>* > recvInterfaces_; //- 3D MUI interfaces to receive from
#endif
        bool receiving_;        //- Determine whether this interface will receive
        DynamicList<Foam::word> recvInterfaceNames_;

        scalar refLength_;
        scalar oneOverRefLength_;
        scalar refVelocity_;
        scalar oneOverRefVelocity_;
        scalar refTime_;
        scalar oneOverRefTime_;
        const List<word>& idList;
        const reducedUnits& rU_;
        List<List<polyMolecule*> > molHistory_;
        List<List<label> > molId_;
#ifdef USE_MUI    	
    	mui::chrono_sampler_exact3d* chrono_sampler;
#endif
    	bool fixedRegion_;
        vector fixedRegionMin_;
        vector fixedRegionMax_;
        bool fixedBounds_;
        vector fixedBoundMin_;
        vector fixedBoundMax_;
        vector fixedBoundNorm_;
        vector fixedBoundZeroThick_;
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

    // Private Member Functions

        polyMolecule* insertMolecule
	(
	    point& position,
	    const label& id,
	    const bool ghost,
	    vector& velocity
	);
        polyMolecule* checkForOverlaps(polyMolecule* newMol, const scalar& potEnergyLimit);

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
        void initialConfiguration();

        //- run at start of time step before molecules move
        void controlAfterMove(label stage);

        //- run after forces are calculated
        void controlAfterForces();

        //- receive molecules that fall into the coupled region
        void receiveCoupledRegion(bool init);

        //- send forces acting on molecules in coupled region
        void sendCoupledRegionForces();

        //- Collate and delete molecules that have crossed a coupling boundary
        void findCoupledMolecules();

        //- send molecules that have crossed a coupling boundary
        void sendCoupledMolecules();

        //- receive molecules that have crossed a coupling boundary
        void receiveCoupledParcels();

        //- Insert parcels that have been received through a coupling boundary
        bool insertCoupledMolecules();

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
        void forget();

        //- Forget the time frame received at specific iteration for all interfaces
        void forget(label iteration);

        //- Forget the time frame received at specific iteration for specific interfaces
        void forget(label iteration, label interface);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
