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
#include "binModel.H"

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

private:

    // Private data

        //- subDictionarys containing the properties
        dictionary propsDict_;
        dictionary propsDictSend_;
        dictionary propsDictRecv_;

        autoPtr<binModel> binModel_;

        List<label> molIds_;

        bool output_;

        //- MUI coupling interfaces (1D & 2D unused in this controller derivation)
        couplingInterface1d oneDInterfaces_;
        couplingInterface2d twoDInterfaces_;
        couplingInterface3d threeDInterfaces_;

#ifdef USE_MUI
        //- Calculated cell centres
        List<mui::point3d> cellCentres_;
#endif
        //- Sending properties
#ifdef USE_MUI
        List<mui::uniface<mui::config_3d>* > sendInterfaces_; //- 3D MUI interfaces to send through
#endif
        bool sending_;          //- Determine whether this controller will send
        List<Foam::word> sendInterfaceNames_;

        //- Receiving properties
#ifdef USE_MUI
        List<mui::uniface<mui::config_3d>* > recvInterfaces_; //- 3D MUI interfaces to receive from
#endif
        bool receiving_;        //- Determine whether this controller will receive
        List<Foam::word> recvInterfaceNames_;

    // Private Member Functions

        template<int cplType>
        void initConfig();

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
            couplingInterface1d &oneDInterfaces,
            couplingInterface2d &twoDInterfaces,
            couplingInterface3d &threeDInterfaces
        );

    // Destructor

        ~dsmcMdCoupling();


    // Member Functions

        //- create an initial configuration
        void initialConfiguration();

        //- run at start of DSMC time step before particles move
        void controlParcelsBeforeMove();

        //- run at end of DSMC time step
        void calculateProperties();

        //- send parcels that fall into the MD coupled region to coupled mdFoamPlus
        void sendCoupledRegion();

        //- send parcels to mdFoamPlus that have crossed a coupling boundary
        void sendCoupledParcels();

        //- receive molecules from mdFoamPlus that have crossed a coupling boundary
        void receiveCoupledMolecules();

        //- output data (optional)
        void output
        (
            const fileName& fixedPathName,
            const List<fileName>& timePaths
        );

        //- update properties from a modified dictionary (optional)
        void updateProperties(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
