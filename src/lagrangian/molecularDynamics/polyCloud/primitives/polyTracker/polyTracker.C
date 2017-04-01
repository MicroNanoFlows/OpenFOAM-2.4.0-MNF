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

#include "polyTracker.H"
#include "polyMoleculeCloud.H"

namespace Foam
{
    

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

polyTracker::polyTracker(polyMoleculeCloud& cloud)
:
    cloud_(cloud),
    id_(-1)
{
    
}







// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyTracker::~polyTracker()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polyTracker::track()
{
    molsOnProc_.clear();
    molTrackingNumbers_.clear();
    
    IDLList<polyMolecule>::iterator mol(cloud_.begin());
    
    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {
        molsOnProc_.append(&mol());
        molTrackingNumbers_.append(mol().trackingNumber());
    }    
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
}

// ************************************************************************* //
