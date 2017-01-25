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
    regionRadii_ = scalarList(propsDict_.lookup("regionRadii"));
    
    nScheduledDestinations_ = destinations_.size();
    
    if( (nScheduledDestinations_ <= 0) || (nScheduledDestinations_ != regionRadii_.size())  )
    {
        FatalError
            << "Something went wrong with the localGoalsSchedule " << endl
            << " lists must match in size (n and n-1) and need to have sizes greater than zero"
            << endl << "destinations = " << destinations_
            << endl << "regionRadii = " << regionRadii_ << endl;
        
    } 
    
    label nextDest = 0;
    
            IDLList<agent>::iterator mol(cloud_.begin());

            for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
            {    
                if(findIndex(agentIds_, mol().id()) != -1)
                {
                    //if(boxes_[nextDest].contains(mol().position())) //WILL THIS WORK FOR AGENTS SEPARATELY?
                    
                        //vector r = destinations_[nextDest] - mol().position(); // NEED AGENT'S POSITION, MAKE SURE IT'S UNIT VECTOR
                        
//                         Info << "Changing destinations of agents of type " << agentIds_ 
//                         << " to destination " << r << endl; 
                        
                        mol().d()= destinations_[nextDest];
                    
                }
            } 
}

// is called every time-step to set new destinations of agents
void localGoalsSchedule::setSchedule()
{
//     timeIndex_ += deltaT_;
    
    // next destination index
//     label nextDest=counter_+1;
    
//     if(nextDest < nScheduledDestinations_)
//     {
//         if(regionRadii_[nextDest] >= timeIndex_)
//         {
            // next destination
            
//             vector r = destinations_[nextDest];
    
            //label index = 0;   // label used for storing int data     
    
            IDLList<agent>::iterator mol(cloud_.begin());

            for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
            {    
                if(findIndex(agentIds_, mol().id()) != -1)
                {
                    
                    label index = 0;
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
                    
                    // check if agent has arrived at its destination
                    if(mag(destinations_[index]-mol().position()) < regionRadii_[index])
                    {
//                         vector r = destinations_[nextDest] - mol().position(); // NEED AGENT'S POSITION, MAKE SURE IT'S UNIT VECTOR
                        
//                         Info << "Changing destinations of agents of type " << agentIds_ 
//                         << " to destination " << r << endl; 
                        if(index < nScheduledDestinations_-1)
                        {
                            mol().d()= destinations_[index+1];
                        }
                    }
                }
            }
            
//             counter_++;        
//     }
    
}



} // End namespace Foam

// ************************************************************************* //
