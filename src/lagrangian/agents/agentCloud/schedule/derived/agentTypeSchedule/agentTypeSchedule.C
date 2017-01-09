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

Description

\*---------------------------------------------------------------------------*/

#include "agentTypeSchedule.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agentTypeSchedule, 0);

addToRunTimeSelectionTable(scheduleModel, agentTypeSchedule, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
agentTypeSchedule::agentTypeSchedule
(
    Time& time,
    agentCloud& cloud,
    const dictionary& dict
)
:
    scheduleModel(time, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    timeIndex_(0.0),
    counter_(-1),
    deltaT_(time_.deltaT().value())
{

    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();    

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentTypeSchedule::~agentTypeSchedule()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void agentTypeSchedule::initialConfiguration()
{
    destinations_ = vectorList(propsDict_.lookup("destinations"));
    scheduledTimes_ = scalarList(propsDict_.lookup("scheduledTimes"));
    
    nScheduledDestinations_ = destinations_.size();
    
    if( (nScheduledDestinations_ <= 0) || (nScheduledDestinations_ != scheduledTimes_.size()) )
    {
        FatalError
            << "Something went wrong with the agentTypeSchedule " << endl
            << " lists must match in size and need to have sizes greater than zero"
            << endl << "destinations = " << destinations_
            << endl << "scheduledTimes = " << scheduledTimes_ << endl;
        
    }    
}

// is called every time-step to set new destinations of agents
void agentTypeSchedule::setSchedule()
{
    timeIndex_ += deltaT_;
    
    // next destination index
    label nextDest=counter_+1;
    
    if(nextDest < nScheduledDestinations_)
    {
        if(scheduledTimes_[nextDest] >= timeIndex_)
        {
            // next destination
            
            vector r = destinations_[nextDest];
            
            Info << "Changing destinations of agents of type " << agentIds_ 
                 << " to destination " << r << endl; 
                 
            IDLList<agent>::iterator mol(cloud_.begin());

            for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
            {
                if(findIndex(agentIds_, mol().id()) != -1)
                {
                    mol().d()=r;
                }
            }
            
            counter_++;        
        }
    }
}



} // End namespace Foam

// ************************************************************************* //
