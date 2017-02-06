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

#include "agentsOnEllipse.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agentsOnEllipse, 0);

addToRunTimeSelectionTable(agentConfiguration, agentsOnEllipse, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
agentsOnEllipse::agentsOnEllipse
(
    agentCloud& molCloud,
    const dictionary& dict
)
:
    agentConfiguration(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))    
{
    treshold_ = 0.001;
    
    checkClosedEndedBorders();
    
    dX_ = 0.5;
    
    if (propsDict_.found("distanceBetweenAgents"))
    {
        dX_ = readScalar(propsDict_.lookup("distanceBetweenAgents"));
    }
    
    R1_ = readScalar(propsDict_.lookup("R1"));
    R2_ = readScalar(propsDict_.lookup("R2"));
    
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentsOnEllipse::~agentsOnEllipse()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void agentsOnEllipse::setInitialConfiguration()
{
//     label initialSize = cloud_.size();

    Info << nl << "Initialising agentsOnEllipse" << endl;

    const word idName(propsDict_.lookup("agentId")); 
    const List<word>& idList(cloud_.cP().agentIds());

    label id = findIndex(idList, idName);

    if(id == -1)
    {
        FatalErrorIn("agentsOnEllipse::setInitialConfiguration()")
            << "Cannot find molecule id: " << idName << nl << "in idList."
            << exit(FatalError);
    }
    
    bool frozen = true;

    if (propsDict_.found("frozen"))
    {
        frozen = Switch(propsDict_.lookup("frozen"));
    }
    
    
   
    DynamicList<vector> agentPositions;


    theta = S/(x(i-1)*x(i-1))
    

    Info << "checking for overlaps... " << endl;
    
    DynamicList<vector> positions;    

    scalar tolerance = 0.05;
    
    forAll(agentPositions, i)
    {
        vector rI = agentPositions[i];
        
        bool overlapping = false;
        
        forAll(positions, j)
        {
            vector rJ = positions[j];
            
            scalar rIJMag = mag(rI - rJ);
            
            if(rIJMag < tolerance)
            {
                overlapping = true;
            }
        }
    
        if (!overlapping)
        {
            positions.append(rI);
        }
    }
 
    
    Info<< "... done. Number of agents to insert = "
        << positions.size() << endl;
        
    if(positions.size() > 1e6)
    {
        FatalErrorIn("agentsOnEllipse::setInitialConfiguration()")
            << "Too many agents to insert"
            << exit(FatalError);
    }
    
    label nAgentsInserted = 0;    
    
    // insert molecules
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
        
        // initialise agent using a random variable 
        vector v = vector::zero;
        
        scalar massI = 1;
        scalar radius = 0.5;
        scalar desiredSpeed = 0.0;
        
        v.z()=0.0;
//         d.add(mag(v));
        
        if(cell != -1)
        {
            insertAgent
            (
                positions[i],
                cell,
                tetFace,
                tetPt,
                id,
                massI,
                radius,
                desiredSpeed,
                frozen,
                v
            );
            
            nAgentsInserted++;
        }
    }
    
    Info << "no of wall-agents inserted = " << nAgentsInserted << endl;
}






} // End namespace Foam

// ************************************************************************* //
