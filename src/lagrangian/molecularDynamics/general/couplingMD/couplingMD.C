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

\*---------------------------------------------------------------------------*/

#include "couplingMD.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::couplingMD::couplingMD(scalar dim, string name)
{
	_appName = name;
	_uri = "mpi://";
	_uri.append(name);
	_uri.append("/ifs");
	_dim = dim;

	#ifdef USE_MUI
		if(_dim == 1)
		{
			_interface1d = new mui::uniface1d(_uri);
		}
		else if (_dim == 2)
		{
			_interface2d = new mui::uniface2d(_uri);
		}
		else if (_dim == 3)
		{
			_interface3d = new mui::uniface3d(_uri);
		}
	#endif
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::couplingMD::~couplingMD()
{
	#ifdef USE_MUI
		if(_dim == 1)
		{
			delete _interface1d;
		}
		else if (_dim == 2)
		{
			delete _interface2d;
		}
		else if (_dim == 3)
		{
			delete _interface3d;
		}
	#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#ifdef USE_MUI
mui::uniface1d* Foam::couplingMD::interface1d() const
{
    return _interface1d;
}

mui::uniface2d* Foam::couplingMD::interface2d() const
{
    return _interface2d;
}

mui::uniface3d* Foam::couplingMD::interface3d() const
{
    return _interface3d;
}
#endif
// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
