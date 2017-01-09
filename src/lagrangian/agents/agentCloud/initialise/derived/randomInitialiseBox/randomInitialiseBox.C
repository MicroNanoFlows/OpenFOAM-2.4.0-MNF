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

#include "randomInitialiseBox.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(randomInitialiseBox, 0);

addToRunTimeSelectionTable(agentConfiguration, randomInitialiseBox, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
randomInitialiseBox::randomInitialiseBox
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

randomInitialiseBox::~randomInitialiseBox()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void randomInitialiseBox::setInitialConfiguration()
{
    label initialSize = cloud_.size();

    Info << nl << "Test initialise " << endl;
    
    
//     const scalar temperature(readScalar(propsDict_.lookup("temperature")));

//     const vector bulkVelocity(propsDict_.lookup("bulkVelocity"));

    const word idName(propsDict_.lookup("agentId")); 
    const List<word>& idList(cloud_.cP().agentIds());

    label id = findIndex(idList, idName);

    if(id == -1)
    {
        FatalErrorIn("randomInitialiseBox::setInitialConfiguration()")
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
    
    label N(readLabel(propsDict_.lookup("N")));
    
    label nAgents = 0;

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
        nAgents++;
    }    

    
    
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

void randomInitialiseBox::setBoundBox
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
