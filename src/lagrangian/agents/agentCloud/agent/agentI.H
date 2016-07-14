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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


inline Foam::agent::agent
(
    const polyMesh& mesh,            
    const agentProperties& cP,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const vector& v,
    const vector& a,
    const vector& f,
    const vector& specialPosition,
    const scalar& R,
    const scalar& frac,   
    const label special,
    const label id,
    const label trackingNumber
)
:
    particle(mesh, position, cellI, tetFaceI, tetPtI),
    v_(v),
    a_(a),
    f_(f),
    specialPosition_(specialPosition),
    R_(R),
    frac_(frac),
    special_(special),
    id_(id),
    trackingNumber_(trackingNumber)
{}

// * * * * * * * * * * * * agent Member Functions * * * * * * * * * * * * //


inline const Foam::vector& Foam::agent::v() const
{
    return v_;
}


inline Foam::vector& Foam::agent::v()
{
    return v_;
}


inline const Foam::vector& Foam::agent::a() const
{
    return a_;
}


inline Foam::vector& Foam::agent::a()
{
    return a_;
}


inline const Foam::vector& Foam::agent::specialPosition() const
{
    return specialPosition_;
}


inline Foam::vector& Foam::agent::specialPosition()
{
    return specialPosition_;
}



inline Foam::label Foam::agent::special() const
{
    return special_;
}

inline Foam::label& Foam::agent::special()
{
    return special_;
}
/*
inline bool Foam::agent::tethered() const
{
    return special_ == SPECIAL_TETHERED;
}*/

inline Foam::label Foam::agent::id() const
{
    return id_;
}

inline Foam::label& Foam::agent::id()
{
    return id_;
}

inline const Foam::scalar& Foam::agent::R() const
{
    return R_;
}

inline Foam::scalar& Foam::agent::R()
{
    return R_;
}

inline const Foam::scalar& Foam::agent::fraction() const
{
    return frac_;
}

inline Foam::scalar& Foam::agent::fraction()
{
    return frac_;
}

inline const Foam::label& Foam::agent::trackingNumber() const
{
    return trackingNumber_;
}

inline Foam::label& Foam::agent::trackingNumber()
{
    return trackingNumber_;
}

inline bool Foam::agent::frozen()
{
    if(special_ == SPECIAL_FROZEN)
    {
        return true;
    }
    else
    {
        return false;
    }
}

inline bool Foam::agent::referred()
{
    return special_ == SPECIAL_REFERRED;
}


// ************************************************************************* //
