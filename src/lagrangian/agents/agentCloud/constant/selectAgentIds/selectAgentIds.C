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
    selectAgentIds

Description

\*----------------------------------------------------------------------------*/

#include "selectAgentIds.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
selectAgentIds::selectAgentIds()
:
    agentIds_()
{}


selectAgentIds::selectAgentIds
(
    const agentProperties& cP,
    const dictionary& dict
)
:
    agentIds_()
{
    const List<word> agents (dict.lookup("agentIds"));
    
    if(agents.size() > 0)
    {
        if
        ( 
            (agents.size() == 1) &&
            (agents[0] == "ALL")
        )
        {
            agentIds_.setSize(cP.nAgentTypes(), -1);
            
            forAll(cP.agentIds(), i)
            {
                agentIds_[i] = i;
            }
        }
        else
        {
            DynamicList<word> agentsReduced(0);
        
            forAll(agents, i)
            {
                const word& agentName(agents[i]);
        
                if(findIndex(agentsReduced, agentName) == -1)
                {
                    agentsReduced.append(agentName);
                }
            }
        
            agentIds_.setSize(agentsReduced.size(), -1);
            agentIdNames_.setSize(agentsReduced.size());
            
            forAll(agentsReduced, i)
            {
                const word& agentName(agentsReduced[i]);
        
                label agentId(findIndex(cP.agentIds(), agentName));
        
                if(agentId == -1)
                {
                    FatalErrorIn
                    (
                        "selectAgentIds::selectAgentIds()"
                    )
                        << "Cannot find id: " << agentName << nl << "in dictionary."
                        << exit(FatalError);
                }
        
                agentIds_[i] = agentId;
                agentIdNames_[i] = agentName;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "selectAgentIds::selectAgentIds()"
        )
            << "agentIds need to be greater than 0 in dictionary."
            << exit(FatalError);
    }
}

selectAgentIds::selectAgentIds
(
    const agentProperties& cP,
    const dictionary& dict,
    const word& agentIdsHeader
)
:
    agentIds_()
{
    const List<word> agents (dict.lookup(agentIdsHeader));

    if(agents.size() > 0)
    {
        if
        ( 
            (agents.size() == 1) &&
            (agents[0] == "ALL")
        )
        {
            agentIds_.setSize(cP.nAgentTypes(), -1);
            
            forAll(cP.agentIds(), i)
            {
                agentIds_[i] = i;
            }
        }
        else
        {
            DynamicList<word> agentsReduced(0);
        
            forAll(agents, i)
            {
                const word& agentName(agents[i]);
        
                if(findIndex(agentsReduced, agentName) == -1)
                {
                    agentsReduced.append(agentName);
                }
            }
        
            agentIds_.setSize(agentsReduced.size(), -1);
            agentIdNames_.setSize(agentsReduced.size());
            
            forAll(agentsReduced, i)
            {
                const word& agentName(agentsReduced[i]);
        
                label agentId(findIndex(cP.agentIds(), agentName));
        
                if(agentId == -1)
                {
                    FatalErrorIn
                    (
                        "selectAgentIds::selectAgentIds()"
                    )
                        << "Cannot find id: " << agentName << nl << "in dictionary."
                        << exit(FatalError);
                }
        
                agentIds_[i] = agentId;
                agentIdNames_[i] = agentName;
                
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "selectAgentIds::selectAgentIds()"
        )
            << "agentIds need to be greater than 0 in dictionary."
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

selectAgentIds::~selectAgentIds()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// void selectAgentIds::operator=(const selectAgentIds& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("selectAgentIds::operator=(const selectAgentIds&)")
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

// Ostream& operator<<(Ostream& os, const selectAgentIds& d)
// {
//     os  << d.binWidth_
//         << static_cast<const Map<label>&>(d);
// 
//     // Check state of Ostream
//     os.check
//     (
//         "Ostream& operator<<(Ostream&, "
//         "const selectAgentIds&)"
//     );
// 
//     return os;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
