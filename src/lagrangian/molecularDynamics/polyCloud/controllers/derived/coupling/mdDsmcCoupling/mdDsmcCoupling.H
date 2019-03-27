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
        scalar refTime_;
        scalar oneOverRefTime_;
        const List<word>& idList;
        const reducedUnits& rU_;
        List<List<polyMolecule*> > molHistory_;
        List<List<bool> > molChanged_;
#ifdef USE_MUI    	
    	mui::chrono_sampler_exact3d* chrono_sampler;
#endif
    	point meshMin_;
		point meshMax_;
		bool refVal_;
		vector refValue_;
		bool fixedBounds_;
		vector fixedBoundMin_;
		vector fixedBoundMax_;
		vector roundCorr_;

    // Private Member Functions

        polyMolecule* insertMolecule
		(
		   const point& position,
		   const label& id,
		   const bool& frozen,
		   const scalar& temperature,
		   vector& bulkVelocity
		);

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

        //- run at start of MD time step before molecules move
        void controlAfterMove(int stage);

        //- run at end of MD time step
        void calculateProperties();

        //- receive molecules that fall into the MD coupled region from coupled dsmcFoamPlus
        void receiveCoupledRegion(bool init);

        //- send molecules to dsmcFoamPlus that have crossed a coupling boundary
        void sendCoupledMolecules();

        //- receive parcels from dsmcFoamPlus that have crossed a coupling boundary
        void receiveCoupledParcels();

        //- output data (optional)
        void output
        (
            const fileName& fixedPathName,
            const List<fileName>& timePaths
        );

        //- update properties from a modified dictionary (optional)
        void updateProperties(const dictionary&);

        //- Create a coupling barrier at current coupling time for all sending interfaces
        void barrier();

        //- Create a coupling barrier at specific coupling time for all sending interfaces
        void barrier(scalar time);

        //- Create a coupling barrier at current coupling time for all sending interfaces
        void barrier(label interface);

        //- Create a coupling barrier at specific coupling time for specific sending interface
        void barrier(scalar time, label interface);

        //- Forget the time frame received at current coupling time for all interfaces
        void forget();

        //- Forget the time frame received at specific time for all interfaces
        void forget(scalar time);

        //- Forget the time frame received at specific time for specific interfaces
        void forget(scalar time, label interface);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
