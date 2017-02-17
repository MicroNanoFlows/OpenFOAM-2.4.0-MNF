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

#include "tiltedCylinder.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tiltedCylinder::tiltedCylinder()
:
    startPoint_(vector::zero),
    endPoint_(vector::zero),
    radius_(0.0)
{
}


Foam::tiltedCylinder::tiltedCylinder
(
    const vector& startPoint,
    const vector& endPoint,
    const scalar& radius    
)
:
    startPoint_(startPoint),
    endPoint_(endPoint),
    radius_(radius)
{
    // test if radius if larger than zero
    
    setBox();
}




// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tiltedCylinder::~tiltedCylinder()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tiltedCylinder::setBox
(
    const vector& startPoint,
    const vector& endPoint,
    const scalar& radius       
)
{
    startPoint_=startPoint;
    endPoint_=endPoint;
    radius_=radius;
    
    setBox();
}

void Foam::tiltedCylinder::setBox()
{
    normal_ = endPoint_ - startPoint_;
    
    magRSE_ = mag(normal_);
    
    normal_ /= magRSE_;
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// void Foam::tiltedCylinder::operator=(const tiltedCylinder& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn
//         (
//             "Foam::tiltedCylinder::operator=(const Foam::tiltedCylinder&)"
//         )   << "Attempted assignment to self"
//             << abort(FatalError);
//     }
// }


// ************************************************************************* //
