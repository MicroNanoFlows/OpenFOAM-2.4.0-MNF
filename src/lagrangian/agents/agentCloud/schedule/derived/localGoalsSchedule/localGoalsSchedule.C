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

#include "localGoalsSchedule.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(localGoalsSchedule, 0);

addToRunTimeSelectionTable(scheduleModel, localGoalsSchedule, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// void localGoalsSchedule::setBoundBoxes()
// {
//     PtrList<entry> boxList(propsDict_.lookup("boxes"));
// 
//     boxes_.setSize(boxList.size());
//     destinations_.setSize(boxList.size());
// 
//     forAll(boxList, b)
//     {
//         const entry& boxI = boxList[b];
//         const dictionary& dict = boxI.dict();
// 
// //         vector startPoint = dict.lookup("startPoint");
// //         vector endPoint = dict.lookup("endPoint");
//         vector goalPoint = dict.lookup("goalPoint");
//         scalar radius = readScalar(
//         destinations_[b] = goalPoint;
//         boxes_[b].resetBoundedBox(startPoint, endPoint);
//     }
// }

// Construct from components
localGoalsSchedule::localGoalsSchedule
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
    
//     setBoundBoxes();

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

localGoalsSchedule::~localGoalsSchedule()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void localGoalsSchedule::initialConfiguration()
{
    destinations_ = vectorList(propsDict_.lookup("destinations"));
    nScheduledDestinations_ = destinations_.size();
    
    regionRadii_ = scalarList(propsDict_.lookup("regionRadii"));
    
    if (propsDict_.found("timeAtDestination"))
    {
        timeAtDestination_ = scalarList(propsDict_.lookup("timeAtDestination"));
    }
    else
    {
        timeAtDestination_.setSize(nScheduledDestinations_, 0.0);
    }

    
    if( (nScheduledDestinations_ <= 0) || 
        (nScheduledDestinations_ != regionRadii_.size()) ||  
        (nScheduledDestinations_ != timeAtDestination_.size())
    )
    {
        FatalError
            << "Something went wrong with the localGoalsSchedule " << endl
            << " lists must match in size and need to have sizes greater than zero"
            << endl << "destinations = " << destinations_
            << endl << "regionRadii = " << regionRadii_
            << endl << "timeAtDestination = " << timeAtDestination_
            << nl << abort(FatalError);  
        
    } 
    
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
void localGoalsSchedule::setSchedule()
{
    IDLList<agent>::iterator mol(cloud_.begin());

    for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
    {    
        if(findIndex(agentIds_, mol().id()) != -1)
        {
            if(mol().t() > 0 )
            {
                mol().t() -= deltaT_;
                
                if(mol().t() < 0)
                {
                    mol().t() = 0.0;
                }
            }
            else
            {
                label index = -1;
                
                // function that produces the index from the agent's d
                forAll(destinations_, i)
                {    
                    if (mag(mol().d() - destinations_[i]) < SMALL) 
                    {
                        index = i;
                    }
                }
                
                if(index == -1)
                {
                    mol().d() = destinations_[0];
                }
                
                // check if agent has arrived at its destination
                if(mag(destinations_[index]-mol().position()) < regionRadii_[index])
                {
                    if(index < nScheduledDestinations_-1)
                    {
                        mol().d()= destinations_[index+1];
                        mol().t()=timeAtDestination_[index];
                    }
                }
            }
        }
    }
}



} // End namespace Foam

// ************************************************************************* //
