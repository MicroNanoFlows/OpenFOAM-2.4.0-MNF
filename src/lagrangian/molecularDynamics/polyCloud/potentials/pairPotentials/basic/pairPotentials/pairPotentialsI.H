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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


inline List< autoPtr<pairPotentialModel> >& pairPotentials::pairPots()
{
    return pairPotentials_;
}

inline const List< autoPtr<pairPotentialModel> >& pairPotentials::pairPots() const
{
    return pairPotentials_;
}

inline const List<word>& pairPotentials::pairNames() const
{
    return pairNames_;
}

inline const List<label>& pairPotentials::pairIds() const
{
    return pairIds_;
}

inline autoPtr<pairPotentialModel>& pairPotentials::electrostaticPot()
{
    return electrostaticPotential_;
}

inline const autoPtr<pairPotentialModel>& pairPotentials::electrostaticPot() const
{
    return electrostaticPotential_;
}

// supply molId I, molId J and i and j pot site indexes (see format of pairPotNames_ in cP_)
//make inline
inline label pairPotentials::getIndexFromPairPotentialSites
(
    const label& idI,
    const label& idJ,
    const label& i,
    const label& j
)
{
    label idsI=cP_.pairPotNames_to_pairPotSitesList()[idI][i];
    label idsJ=cP_.pairPotNames_to_pairPotSitesList()[idJ][j];
    return pairPotIdList_to_pairPotentials_[idsI][idsJ];
}

// supply molId I, molId J and sI and sJ site indexes (see format of siteNames_ in cP_)
//make inline

inline label pairPotentials::getIndexFromSites
(
    const label& idI,
    const label& idJ,
    const label& sI,
    const label& sJ
)
{
    label idsI=cP_.siteNames_to_siteIdList()[idI][sI];
    label idsJ=cP_.siteNames_to_siteIdList()[idJ][sJ];
    return siteIdList_to_pairPotentials_[idsI][idsJ];
}

// inline autoPtr<pairPotentialModel>& pairPotentials::pairPot
// (
//     const label pairSiteA,
//     const label pairSiteB
// )
// {
//      label k = pairPotLinks_[pairSiteA][pairSiteB];
//     
//      return pairPotentials_[k];
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
