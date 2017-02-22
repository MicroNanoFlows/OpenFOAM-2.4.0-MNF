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

#include "agentNumbers.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agentNumbers, 0);

addToRunTimeSelectionTable(agentController, agentNumbers, dictionary);


void agentNumbers::setBoundBox
(
    const dictionary& propsDict,
    boundedBox& bb,
    const word& name 
)
{
    const dictionary& dict(propsDict.subDict(name));
    
    vector startPoint = dict.lookup("startPoint");
    vector endPoint = dict.lookup("endPoint");

    bb.resetBoundedBox(startPoint, endPoint);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
agentNumbers::agentNumbers
(
    Time& t,
    agentCloud& cloud,
    const dictionary& dict
)
:
    agentController(t,  cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    agentId_()
{
    writeInTimeDir_ = false;
    writeInCase_ = false;

    {
        const List<word>& idList(cloud_.cP().agentIds());
        const word agentId = propsDict_.lookup("agentId");
        agentId_ = findIndex(idList, agentId);
    
        if(agentId_ == -1)
        {
            FatalErrorIn("agentNumbers::agentNumbers()")
                << "Cannot find agentId: " << agentId << nl << "in: "
                << time_.time().system()/"controllersDict"
                << exit(FatalError);
        }
    }
    
    N_ = readLabel(propsDict_.lookup("N"));
    
    setBoundBox(propsDict_, samplingBox_, "samplingBox");
    setBoundBox(propsDict_, controlBox_, "controlBox");
    
    mass_ = 70;
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentNumbers::~agentNumbers()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void agentNumbers::initialConfiguration()
{}

void agentNumbers::controlBeforeVelocityI()
{}

void agentNumbers::controlBeforeMove()
{}

void agentNumbers::controlBeforeForces()
{}

void agentNumbers::controlDuringForces
(
    agent* molI,
    agent* molJ
)
{}

void agentNumbers::controlAfterForces()
{
    Info << "agentNumbers: control" << endl;

    label nAgents = 0;
    

    IDLList<agent>::iterator p(cloud_.begin());

    for (p = cloud_.begin(); p != cloud_.end(); ++p)
    {
        if(agentId_ ==  p().id())
        {
            if(samplingBox_.contains(p().position()))
            {
                nAgents++;
            }
        }
    }
    
/*    if(Pstream::parRun())
    {
        reduce(nAgents, sumOp<label>());
    }*/        
        
    // number of agents to insert/delete
    
    label N = N_ - nAgents;
    
    if(N > 0)
    {
        Info << "number of agents to insert = " << N  << endl;
        
        insertAgents(N);
    }
}

void agentNumbers::insertAgents(const label& N)
{
    boundedBox& bb = controlBox_;
    
    DynamicList<vector> positions;
    
    for (label i = 0; i < N; i++)
    {
        vector pos = vector
                     ( 
                        (bb.max().x() - bb.min().x())*cloud_.rndGen().scalar01() + bb.min().x(), 
                        (bb.max().y() - bb.min().y())*cloud_.rndGen().scalar01() + bb.min().y(),  
                        0.0
                     );
        
        positions.append(pos);
    }    

    //*** deal with parallel processing here ***//
    
    
    
    label nAgentsInserted = 0;
    
    // insert agents
    
    forAll(positions, i)
    {
        label cell = -1;
        label tetFace = -1;
        label tetPt = -1;

        mesh_.findCellFacePt
        (
            positions[i],
            cell,
            tetFace,
            tetPt
        );
        
        if(cell != -1)
        {
            cloud_.createAgent
            (
                positions[i],
                cell,
                tetFace,
                tetPt,     
                vector::zero, //v
                vector::zero,
                vector::zero,
                vector::zero,
                mass_, // mass 
                0.35,
                1.0, 
                0.0,
                GREAT,
                1.0,
                0.0,
                agentId_,
                0,
                cloud_.getTrackingNumber()
            );
            
            nAgentsInserted++;
        }
    }
    
    Info << "inserted = " << nAgentsInserted << endl;
}


void agentNumbers::controlAfterVelocityII()
{}

void agentNumbers::calculateProperties()
{}

void agentNumbers::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
    }
}

// void agentNumbers::updateProperties(const dictionary& newDict)
// {
//     //- the main controller properties should be updated first
//     updateStateControllerProperties(newDict);
// 
//     propsDict_ = newDict.subDict(typeName + "Properties");
// 
//     model_->updateProperties(propsDict_);
// 
//     readProperties();
// }

/*
void agentNumbers::readProperties()
{

}*/



} // End namespace Foam

// ************************************************************************* //
