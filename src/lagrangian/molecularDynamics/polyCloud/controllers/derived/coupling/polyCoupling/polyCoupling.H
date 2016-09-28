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
    polyCoupling

Description


SourceFiles
    polyCoupling.C

\*---------------------------------------------------------------------------*/

#ifndef polyCoupling_H
#define polyCoupling_H

#include "polyCouplingController.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "binModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyCoupling Declaration
\*---------------------------------------------------------------------------*/

class polyCoupling
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

        bool output_;

        //- MUI coupling interfaces (1D & 2D unused in this controller derivation)
        List<couplingInterface1d> oneDCouplings_;
        List<couplingInterface2d> twoDCouplings_;
		List<couplingInterface3d> threeDCouplings_;

#ifdef USE_MUI
		//- Pointer to correct MUI interface
		mui::uniface3d *currInterface;
		//- Calculated cell centres
		List<mui::point3d> cellCentres_;
#endif

		word cplInterfaceName_; //Name for interface (set in system/couplingDict)
		bool sending_;			//Determine whether this interface will send
		bool receiving_;		//Determine whether this interface will receive

		//- Sending properties
		bool sendMass_;
		bool sendDensity_;

		//- Receiving properties
		List<word> recvInterfaces_;
		bool recvMass_;
		bool recvDensity_;
		List<List<scalar> > recvMassValues_;
		List<List<scalar> > recvDensityValues_;

    // Private Member Functions

		template<int cplType>
		void initConfig();

public:

    //- Runtime type information
    TypeName("polyCoupling");

    // Constructors

        //- Construct from components
    	polyCoupling
        (
            Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict,
			List<couplingInterface1d>& oneDCouplings,
			List<couplingInterface2d>& twoDCouplings,
			List<couplingInterface3d>& threeDCouplings
        );

    // Destructor

        ~polyCoupling();


    // Member Functions

        //- create an initial configuration
        void initialConfiguration();
        
        //- send any coupled data to the MUI interface (end of time-step)
        void sendCoupling();

        //- receive any coupled data to the MUI interface (end of time-step after sendCoupling)
		void receiveCoupling();
        
        //- output data (optional)
        void output
        (
            const fileName& fixedPathName, 
            const fileName& timePath
        );
        
        //- update properties from a modified dictionary (optional)
        void updateProperties(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
