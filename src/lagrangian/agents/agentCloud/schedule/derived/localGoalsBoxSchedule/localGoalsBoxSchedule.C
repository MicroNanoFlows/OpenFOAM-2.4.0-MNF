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

#include "localGoalsBoxSchedule.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(localGoalsBoxSchedule, 0);

addToRunTimeSelectionTable(scheduleModel, localGoalsBoxSchedule, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void localGoalsBoxSchedule::setBoundBoxes()
{
    PtrList<entry> boxList(propsDict_.lookup("boxes"));

    boxes_.setSize(boxList.size());
    destinations_.setSize(boxList.size());
    nScheduledDestinations_ = destinations_.size();
    
    if (propsDict_.found("timeAtDestination"))
    {
        timeAtDestination_ = scalarList(propsDict_.lookup("timeAtDestination"));
    }
    else
    {
        timeAtDestination_.setSize(nScheduledDestinations_, 0.0);
    }    
    
    forAll(boxList, b)
    {
        const entry& boxI = boxList[b];
        const dictionary& dict = boxI.dict();

        vector startPoint = dict.lookup("startPoint");
        vector endPoint = dict.lookup("endPoint");
        boxes_[b].resetBoundedBox(startPoint, endPoint);
        destinations_[b] = boxes_[b].midpoint();        
    }
    
    if( (nScheduledDestinations_ <= 0) || 
        (nScheduledDestinations_ != timeAtDestination_.size())
    )
    {
        FatalError
            << "Something went wrong with the localGoalsSchedule " << endl
            << " lists must match in size and need to have sizes greater than zero"
            << endl << "destinations = " << destinations_
            << endl << "timeAtDestination = " << timeAtDestination_
            << nl << abort(FatalError);  
        
    }     
}

// Construct from components
localGoalsBoxSchedule::localGoalsBoxSchedule
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

localGoalsBoxSchedule::~localGoalsBoxSchedule()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void localGoalsBoxSchedule::initialConfiguration()
{
    setBoundBoxes();    
    
    label nextDest = 0;
    
    IDLList<agent>::iterator mol(cloud_.begin());

    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {    
        if(findIndex(agentIds_, mol().id()) != -1)
        {
            mol().d()= destinations_[nextDest];
        }
    } 
}

// is called every time-step to set new destinations of agents
void localGoalsBoxSchedule::setSchedule()
{
    
    IDLList<agent>::iterator mol(cloud_.begin());

    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {    
        if(findIndex(agentIds_, mol().id()) != -1)
        {
            if(mol().t() > 0 )
            {
		mol().v() = vector::zero;
                mol().t() -= deltaT_;
                
                if(mol().t() < 0)
                {
                    mol().t() = 0.0;
                }
            }
            else
            {
            
                label index = -1;
                //if(boxes_[nextDest].contains(mol().position())) //WILL THIS WORK FOR AGENTS SEPARATELY?
                
                // function that produces the index from the agent's d
                //for (label i = 0; i < nScheduledDestinations_; ++i) //CORRECT????
                forAll(destinations_, i)
                {    
                    if (mag(mol().d() - destinations_[i]) < SMALL) //How to get magnitude of vectors?
                    {
                        index = i;
                    }
                }
                
                if(index == -1)
                {
                    mol().d() = destinations_[0];
                }            
                
                // check if agent has arrived at its destination
                if(boxes_[index].contains(mol().position()))
                {
                    if(index < nScheduledDestinations_-1)
                    {
                        mol().d() = destinations_[index+1];
                        mol().t()=timeAtDestination_[index];
                    }
                }
            }
        }
    }

}



} // End namespace Foam

// ************************************************************************* //
