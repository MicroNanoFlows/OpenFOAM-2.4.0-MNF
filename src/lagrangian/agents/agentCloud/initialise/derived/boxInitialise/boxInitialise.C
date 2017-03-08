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

#include "boxInitialise.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(boxInitialise, 0);

addToRunTimeSelectionTable(agentConfiguration, boxInitialise, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
boxInitialise::boxInitialise
(
    agentCloud& molCloud,
    const dictionary& dict
)
:
    agentConfiguration(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))    
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

boxInitialise::~boxInitialise()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void boxInitialise::setInitialConfiguration()
{
//     label initialSize = cloud_.size();

    Info << nl << "Box initialise " << endl;
    
    
//     const scalar temperature(readScalar(propsDict_.lookup("temperature")));

//     const vector bulkVelocity(propsDict_.lookup("bulkVelocity"));

    const word idName(propsDict_.lookup("agentId")); 
    const List<word>& idList(cloud_.cP().agentIds());

    label id = findIndex(idList, idName);

    if(id == -1)
    {
        FatalErrorIn("boxInitialise::setInitialConfiguration()")
            << "Cannot find molecule id: " << idName << nl << "in idList."
            << exit(FatalError);
    }

    scalar meanMass = readScalar(propsDict_.lookup("meanMass"));
    scalar massRange = readScalar(propsDict_.lookup("massRange"));

    scalar meanRadius = readScalar(propsDict_.lookup("meanRadius"));
    scalar radiusRange = readScalar(propsDict_.lookup("radiusRange"));

    scalar meanDesSpeed = readScalar(propsDict_.lookup("meanDesiredSpeed"));
    scalar desSpeedRange = readScalar(propsDict_.lookup("desiredSpeedRange"));
   
    scalar minV = readScalar(propsDict_.lookup("minV"));
    scalar maxV = readScalar(propsDict_.lookup("maxV"));        
    
    bool frozen = false;

    if (propsDict_.found("frozen"))
    {
        frozen = Switch(propsDict_.lookup("frozen"));
    }

    const scalar binWidth(readScalar(propsDict_.lookup("binWidth")));
    
    distribution d(binWidth);
   
    
    // Bounding box 
    boundedBox bb;
    
    setBoundBox(propsDict_, bb, "boundBox");
    
    label Nx(readLabel(propsDict_.lookup("Nx")));
    label Ny(readLabel(propsDict_.lookup("Ny")));
//     label Nx = 10;
//     label Ny = 10;
    
  
    scalar Lx = bb.span().x();
    scalar Ly = bb.span().y();
    
    scalar DX = Lx / Nx; // spacing
    scalar DY = Ly / Ny;
    
    label nAgents = 0;

    DynamicList<vector> positions;
    
    for (label i = 0; i < Nx; i++)
    {
        for (label j = 0; j < Ny; j++)
        {
            vector pos = bb.min() + vector(1, 0, 0)*i*DX +  vector(1, 0, 0)*DX*0.5 + 
                                    vector(0, 1, 0)*j*DY +  vector(0, 1, 0)*DY*0.5;
                                    
            if(bb.contains(pos))
            {
                positions.append(pos);
                nAgents++;
            }
        }
    }    
    
//     label N(readLabel(propsDict_.lookup("N")));


//     Info << "positions = " << positions << endl;
    
    
    Info << "number of agents, N = " << positions.size() << endl;
    
    label nMolsInserted = 0;
    
    
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
        vector v = setInitialVelocity(minV, maxV);
        
        scalar massI = gaussianDistribution(meanMass, massRange);
        scalar radius = gaussianDistribution(meanRadius, radiusRange);
        scalar desiredSpeed = gaussianDistribution(meanDesSpeed, desSpeedRange);
        
        d.add(massI);
        
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
            
            nMolsInserted++;
        }
    }
    
    
    List< Pair<scalar> > histogram = d.raw();
    
    {    
        OFstream file("distribution.xy");

        if(file.good())
        {
            forAll(histogram, i)
            {
                file 
                    << histogram[i].first() << "\t"
                    << histogram[i].second() << "\t"
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void pairPotentialModel::write()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
   
}

void boxInitialise::setBoundBox
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




} // End namespace Foam

// ************************************************************************* //
