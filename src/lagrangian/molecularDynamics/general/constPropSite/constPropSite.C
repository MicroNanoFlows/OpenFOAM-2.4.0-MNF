/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    constPropSite

Description

\*----------------------------------------------------------------------------*/

#include "constPropSite.H"
#include "graph.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
// void constPropSite::setRadius()
// {
//     for(label i = 0; i < noOfBins_; i++)
//     {
//        radius_[i] = (0.5 + scalar(i)) * binWidth();
//     }
// }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
constPropSite::constPropSite()
:
    siteReferencePosition_(vector::zero),
    siteMass_(0.0),
    siteCharge_(0.0),
    siteId_(0),
    name_(),
    pairPotentialSite_(false),
    electrostaticSite_(false)
{}



constPropSite::constPropSite
(
    const vector& siteReferencePosition,
    const scalar& siteMass,
    const scalar& siteCharge,
    const label& siteId,
    const word& name,
    const bool& pairPotentialSite,
    const bool& electrostaticSite
)
:
    siteReferencePosition_(siteReferencePosition),
    siteMass_(siteMass),
    siteCharge_(siteCharge),
    siteId_(siteId),
    name_(name),
    pairPotentialSite_(pairPotentialSite),
    electrostaticSite_(electrostaticSite)
{
  
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constPropSite::~constPropSite()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// void constPropSite::operator=(const constPropSite& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("constPropSite::operator=(const constPropSite&)")
//             << "Attempted assignment to self"
//             << abort(FatalError);
//     }
// 
//     Map<label>::operator=(rhs);
// 
//     binWidth_ = rhs.binWidth();
// }


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// Ostream& operator<<(Ostream& os, const constPropSite& d)
// {
//     os  << d.binWidth_
//         << static_cast<const Map<label>&>(d);
// 
//     // Check state of Ostream
//     os.check
//     (
//         "Ostream& operator<<(Ostream&, "
//         "const constPropSite&)"
//     );
// 
//     return os;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
