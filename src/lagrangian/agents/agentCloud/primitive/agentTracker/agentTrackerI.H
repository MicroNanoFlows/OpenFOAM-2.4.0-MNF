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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{
    
    
    

inline bool agentTracker::isTrackingNumberAvailable
(
    const label& trackingNumber
)
{
    
    id_ = findIndex(agentTrackingNumbers_, trackingNumber);
    
    if(id_ != -1)
    {
        return true;
    }
    else
    {
        return false;
    }
    
    // deal with parallel processing here
}

inline agent* agentTracker::getAgent()
{
    agent* p = NULL;
    
    if(id_ != -1)
    {
       p = agentsOnProc_[id_];
    }
    
    id_ = -1;
    
    return p;
}

inline const DynamicList<label>& agentTracker::agentTrackingNumbers() const 
{
    return agentTrackingNumbers_;
}

inline DynamicList<agent*>& agentTracker::agentsOnProc()
{
    return agentsOnProc_;
}

inline const DynamicList<agent*>& agentTracker::agentsOnProc() const
{
    return agentsOnProc_;
}

inline label agentTracker::noOfAgentsOnProc()
{
    return agentsOnProc_.size();
}


// ************************************************************************* //
}