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

#include "lattice.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(lattice, 0);

addToRunTimeSelectionTable(configuration, lattice, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
lattice::lattice
(
    molecules& cloud,
    const dictionary& dict
)
:
    configuration(cloud, dict)
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

lattice::~lattice()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void lattice::spawn()
{
    Info << nl << "Spawning: " << type() << endl;
    
    const List<word> types(dict_.lookup("types")); 

    forAll(types, i)
    {
        if(findIndex(cloud_.typeNames(), types[i]) == -1)
        {
            FatalErrorIn("lattice::spawn()")
                << "Cannot find atom type: " << types[i] 
                << nl << "in types list."
                << exit(FatalError);            
        }
    }

    List<vector> sitePositions;
    
    label nSites=0;

    const word latticeType(dict_.lookup("latticeType"));

    if (latticeType=="SC")
    {
        nSites=1;
        sitePositions.setSize(nSites,vector::zero);
        sitePositions[0]=vector(0.0,0.0,0.0);
    }
    else if (latticeType=="BCC")
    {
        nSites=2;
        sitePositions.setSize(nSites,vector::zero);
        sitePositions[0]=vector(0.0,0.0,0.0);
        sitePositions[1]=vector(0.5,0.5,0.5);
    }
    else if (latticeType=="FCC")
    {
        nSites=4;
        sitePositions.setSize(nSites,vector::zero);
        sitePositions[0]=vector(0.0,0.0,0.0);
        sitePositions[1]=vector(0.5,0.5,0.0);
        sitePositions[2]=vector(0.5,0.0,0.5);
        sitePositions[3]=vector(0.0,0.5,0.5);
    }
//     else if (latticeType=="HPC")
//     {
//         nSites=4;
//         sitePositions.setSize(nSites,vector::zero);
//         sitePositions[0]=vector(0.0,0.0,0.0);
//         sitePositions[1]=vector(0.5,0.5,0.0);
//         sitePositions[2]=vector(0.5,0.0,0.5);
//         sitePositions[3]=vector(0.0,0.5,0.5);
//     }
    else
    {
        FatalErrorIn("lattice::spawn()")
            << "The lattice option: " << latticeType 
            << ", is not supported"
            << exit(FatalError);
    }
    
    // Bounding box 

    boundsBox bb;
    
    setBoundsBox(dict_, bb, "boundsBox");

    scalar s(readScalar(dict_.lookup("unitCellSize")));

    label nX = (bb.span().x()/s) + 1;
    label nY = (bb.span().y()/s) + 1;
    label nZ = (bb.span().z()/s) + 1;

    label nAtoms= 0;
    DynamicList<vector> positions;
    
    for (label k = 0; k < nX; k++)
    {
        for (label j = 0; j < nY; j++)
        {
            for (label i = 0; i < nZ; i++)
            {
                for (label iS = 0; iS < nSites; iS++)
                {
                    vector sP = sitePositions[iS];
                    
                    vector pos = vector(1, 0, 0)*(k+sP.x())*s + 
                                 vector(0, 1, 0)*(j+sP.y())*s + 
                                 vector(0, 0, 1)*(i+sP.z())*s + 
                                 bb.min();
                    
                    if(bb.contains(pos))
                    {
                        positions.append(pos);
                        nAtoms++;
                    }
                }
            }
        }
    }
    
    positions.shrink();

//     Info << nl << " No of sites found = " << positions.size() << endl;

    // check for overlaps here (to add)
    
    
    
    forAll(positions, i)
    {
        cloud_.positions().append(positions[i]);
        
        label id = cloud_.getType(types[0]);
        
        if(id != -1)
        {
            cloud_.types().append(id);
        }
        else
        {
            // error
        }
        
        cloud_.charges().append(0.0);
    }
    
}

void lattice::setBoundsBox
(
    const dictionary& propsDict,
    boundsBox& bb,
    const word& name 
)
{
    const dictionary& dict(propsDict.subDict(name));
    
    vector startPoint = dict.lookup("startPoint");
    vector endPoint = dict.lookup("endPoint");

    bb.resetBoundedBox(startPoint, endPoint);
}




} // End namespace Foam

// ************************************************************************* //
