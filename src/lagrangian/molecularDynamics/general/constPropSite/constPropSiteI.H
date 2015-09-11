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
    distribution

Description

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const vector& constPropSite::siteReferencePosition() const
{
    return siteReferencePosition_;
}

inline vector& constPropSite::siteReferencePosition()
{
    return siteReferencePosition_;
}

inline const scalar& constPropSite::siteMass() const
{
    return siteMass_;
}

inline scalar& constPropSite::siteMass()
{
    return siteMass_;
}

inline const scalar& constPropSite::siteCharge() const
{
    return siteCharge_;
}

inline scalar& constPropSite::siteCharge()
{
    return siteCharge_;
}



inline const label& constPropSite::siteId() const
{
    return siteId_;
}


inline label& constPropSite::siteId()
{
    return siteId_;
}

inline const word& constPropSite::name() const
{
    return name_;
}

inline word& constPropSite::name()
{
    return name_;
}

inline const bool& constPropSite::pairPotentialSite() const
{
    return pairPotentialSite_;
}

inline bool& constPropSite::pairPotentialSite()
{
    return pairPotentialSite_;
}


inline const bool& constPropSite::electrostaticSite() const
{
    return electrostaticSite_;
}

inline bool& constPropSite::electrostaticSite()
{
    return electrostaticSite_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
