/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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
    Foam::couplingMD

Description

SourceFiles
    couplingMD.C

\*---------------------------------------------------------------------------*/

#ifndef couplingMD_H
#define couplingMD_H

#include "scalar.H"

// NOTE:
// MUI header included if the switch -DUSE_MUI included during compilation.
#ifdef USE_MUI
	#include "mui.h"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coupling Declaration
\*---------------------------------------------------------------------------*/

class couplingMD
{
	// Private data
private:
	string _appName;
	string _uri;
	scalar _dim;

#ifdef USE_MUI
	mui::uniface1d *_interface1d;
	mui::uniface2d *_interface2d;
	mui::uniface3d *_interface3d;
#endif

public:

    // Constructors

        //- Construct
		couplingMD
		(
			scalar dimensions,
			string name
		);

    // Destructor

        ~couplingMD();

    // Member Functions
#ifdef USE_MUI
        mui::uniface1d* interface1d() const;
        mui::uniface2d* interface2d() const;
        mui::uniface3d* interface3d() const;
#endif
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct couplingInterface
{
	word couplingName;
	word interfaceName;
	couplingMD* interface;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
