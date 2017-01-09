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

inline bool agentMicroElements::run()
{
    bool end = true;
    
    if(macroTime_ >= macroEndTime_)
    {
        end = false;
    }
    
    return end;
}

inline const label& agentMicroElements::N() const
{
    return N_;
}
    

inline List< autoPtr<agentCloud> >& agentMicroElements::clouds()
{
	return clouds_;
}

inline const List< autoPtr<agentCloud> >& agentMicroElements::clouds() const
{
	return clouds_;
}

inline label agentMicroElements::iMinus (const label& i, const label& value)
{
//     label k = i-1;
    
    label k = i-value;
    
    if(k < 0)
    {
        k = nMicro_ + k;
    }
    
    return k;
}

inline label agentMicroElements::iPlus (const label& i, const label& value)
{
//     label k = i+1;
    
    label k = i+value;
    
    if(k >= nMicro_)
    {
        k = k - nMicro_;
    }
    
    return k;
}



inline const scalar& agentMicroElements::macroTime() const
{
    return macroTime_;
}

inline const scalar& agentMicroElements::microTime(const label& i) const
{
    return microTimes_[i];
}



// inline const bool& agentMicroElements::pseudo() const
// {
//     return pseudo_;
// }

inline const labelListList& agentMicroElements::procsMDaddressing() const
{
    return procsMDaddressing_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
