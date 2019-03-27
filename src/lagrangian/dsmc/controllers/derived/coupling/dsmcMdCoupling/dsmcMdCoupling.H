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
    dsmcMdCoupling

Description


SourceFiles
    dsmcMdCoupling.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcMdCoupling_H
#define dsmcMdCoupling_H

#include "dsmcCouplingController.H"
#include "vector.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcMdCoupling Declaration
\*---------------------------------------------------------------------------*/

class dsmcMdCoupling
:
    public dsmcCouplingController
{

private:

    // Private data

        //- subDictionarys containing the properties
        dictionary propsDict_;
        dictionary propsDictSend_;
        dictionary propsDictRecv_;

        List<word> typeNames_;
        bool output_;

#ifdef USE_MUI
        //- Calculated cell centres
        List<mui::point3d> cellCentres_;
#endif
        //- Sending properties
#ifdef USE_MUI
        DynamicList<mui::uniface<mui::config_3d>* > sendInterfaces_; //- 3D MUI interfaces to send through
#endif
        bool sending_;          //- Determine whether this controller will send
        DynamicList<word> sendInterfaceNames_;

        //- Receiving properties
#ifdef USE_MUI
        DynamicList<mui::uniface<mui::config_3d>* > recvInterfaces_; //- 3D MUI interfaces to receive from
#endif
        bool receiving_;        //- Determine whether this controller will receive
        DynamicList<word> recvInterfaceNames_;

        List<List<label> > parcelsInCellHistory_;
        scalar refLength_;
		scalar oneOverRefLength_;
		scalar refTime_;
		scalar oneOverRefTime_;
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

    // Private Member Functions

        void insertParcel
        (
            const point& position,
        	const vector& U,
        	const label& typeId
        );

public:

    //- Runtime type information
    TypeName("dsmcMdCoupling");

    // Constructors

        //- Construct from components
        dsmcMdCoupling
        (
            Time& t,
            dsmcCloud& cloud,
            const dictionary& dict,
            couplingInterface2d& twoDInterfaces,
            couplingInterface3d& threeDInterfaces
        );

    // Destructor

        ~dsmcMdCoupling();


    // Member Functions

        //- create an initial configuration
        void initialConfiguration();

        //- run at end of DSMC time step
        void calculateProperties();

        //- send parcels that fall into the MD coupled region to coupled mdFoamPlus
        void sendCoupledRegion(bool init);

        //- send parcels to mdFoamPlus that have crossed a coupling boundary
        void sendCoupledParcels();

        //- receive molecules from mdFoamPlus that have crossed a coupling boundary
        bool receiveCoupledMolecules();

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
