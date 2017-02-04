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

#include "agentsOnBorders.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agentsOnBorders, 0);

addToRunTimeSelectionTable(agentConfiguration, agentsOnBorders, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
agentsOnBorders::agentsOnBorders
(
    agentCloud& molCloud,
    const dictionary& dict
)
:
    agentConfiguration(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))    
{
    borderList_ = List<vectorList>(propsDict_.lookup("bordersList"));

    treshold_ = 0.001;
    
    checkClosedEndedBorders();
    
    dX_ = 0.5;
    
    if (propsDict_.found("distanceBetweenAgents"))
    {
        dX_ = readScalar(propsDict_.lookup("distanceBetweenAgents"));
    }     
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentsOnBorders::~agentsOnBorders()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void agentsOnBorders::checkClosedEndedBorders()
{
    forAll(borderList_, i)
    {
        const vector& rI = borderList_[i][0];
        const vector& rJ = borderList_[i][borderList_[i].size()-1];
        scalar rD = mag(rI - rJ);
        
        if(rD > treshold_)
        {
            FatalErrorIn("reflectiveRectangularBorder::checkClosedEndedBorders()") 
                << "    This border has its end points not connected: " << nl
                << borderList_[i] << nl
                << "    Make the start point and end points the same." 
                << nl << abort(FatalError);  
        }
    }
}

void agentsOnBorders::setInitialConfiguration()
{
//     label initialSize = cloud_.size();

    Info << nl << "Initialising agentsOnBorders" << endl;

    const word idName(propsDict_.lookup("agentId")); 
    const List<word>& idList(cloud_.cP().agentIds());

    label id = findIndex(idList, idName);

    if(id == -1)
    {
        FatalErrorIn("agentsOnBorders::setInitialConfiguration()")
            << "Cannot find molecule id: " << idName << nl << "in idList."
            << exit(FatalError);
    }
    
    bool frozen = true;

    if (propsDict_.found("frozen"))
    {
        frozen = Switch(propsDict_.lookup("frozen"));
    }
    
    
   
    DynamicList<vector> agentPositions;
    
    forAll(borderList_, i)
    {
        label nEdges = borderList_[i].size();
        
        for (int p = 0; p < nEdges-1; p++)
        {
            const vector& v1=borderList_[i][p];
            const vector& v2=borderList_[i][p+1];
            
            vector n = v2-v1;
            scalar mag21 = mag(n);
            n /= mag21;
            
            label nPts = label(mag21/dX_);
            
            // new DX to include round off errors
            scalar dX = mag21/nPts;
            
//             Info << "v1 = " << v1
//             << ", v2 = " << v2
//             << ", mag = " << mag21
//             << " nPts = " << nPts << endl;
            
            for (int j = 0; j < nPts - 1; j++)
            {
                vector point  = v1 + dX*n*j;
                
                agentPositions.append(point);
            }
        }
    }

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
        FatalErrorIn("agentsOnBorders::setInitialConfiguration()")
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
