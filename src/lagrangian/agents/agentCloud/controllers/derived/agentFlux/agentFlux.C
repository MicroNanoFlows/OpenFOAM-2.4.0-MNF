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

#include "agentFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agentFlux, 0);

addToRunTimeSelectionTable(agentController, agentFlux, dictionary);


void agentFlux::setBoundBox
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
agentFlux::agentFlux
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
            FatalErrorIn("agentFlux::agentFlux()")
                << "Cannot find agentId: " << agentId << nl << "in: "
                << time_.time().system()/"controllersDict"
                << exit(FatalError);
        }
    }
    
    N_ = readLabel(propsDict_.lookup("noOfAgents_everyWriteInterval"));

    meanMass_ = readScalar(propsDict_.lookup("meanMass"));
    massRange_ = readScalar(propsDict_.lookup("massRange"));

    meanRadius_ = readScalar(propsDict_.lookup("meanRadius"));
    radiusRange_ = readScalar(propsDict_.lookup("radiusRange"));

    meanDesSpeed_ = readScalar(propsDict_.lookup("meanDesiredSpeed"));
    desSpeedRange_ = readScalar(propsDict_.lookup("desiredSpeedRange"));
    
//     setBoundBox(propsDict_, samplingBox_, "samplingBox");
    setBoundBox(propsDict_, controlBox_, "controlBox");
    
    // better way of adding mass 
//     mass_ = 70;
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentFlux::~agentFlux()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void agentFlux::initialConfiguration()
{}

void agentFlux::controlBeforeVelocityI()
{}

void agentFlux::controlBeforeMove()
{}

void agentFlux::controlBeforeForces()
{}

void agentFlux::controlDuringForces
(
    agent* molI,
    agent* molJ
)
{}

void agentFlux::controlAfterForces()
{
    if(time_.outputTime())
    {
        Info << "agentFlux: control" << endl;
        
        insertAgents(N_);
    }
}

void agentFlux::insertAgents(const label& N)
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
    
    
    scalar massI = gaussianDistribution(meanMass_, massRange_);
    scalar radius = gaussianDistribution(meanRadius_, radiusRange_);
    scalar desiredSpeed = gaussianDistribution(meanDesSpeed_, desSpeedRange_);
    
    
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
                massI, // mass 
                radius,
                desiredSpeed, 
                0.0,
                GREAT,
                1.0,
                0.0,
                agentId_,
                cloud_.getTrackingNumber()
            );
            
            nAgentsInserted++;
        }
    }
    
    Info << "inserted = " << nAgentsInserted << endl;
}


void agentFlux::controlAfterVelocityII()
{}

void agentFlux::calculateProperties()
{}

void agentFlux::output
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

scalar agentFlux::gaussianDistribution
(
    const scalar& mean,
    const scalar& range
 
)
{
    return mean + 0.25*range*cloud_.rndGen().GaussNormal();
}


} // End namespace Foam

// ************************************************************************* //
