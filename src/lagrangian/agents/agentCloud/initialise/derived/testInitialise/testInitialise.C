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

#include "testInitialise.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(testInitialise, 0);

addToRunTimeSelectionTable(agentConfiguration, testInitialise, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
testInitialise::testInitialise
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

testInitialise::~testInitialise()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void testInitialise::setInitialConfiguration()
{
    label initialSize = cloud_.size();

    Info << nl << "Test initialise " << endl;
    
    
    const scalar temperature(readScalar(propsDict_.lookup("temperature")));

//     const vector bulkVelocity(propsDict_.lookup("bulkVelocity"));

    const word idName(propsDict_.lookup("id")); 
    const List<word>& idList(cloud_.cP().agentIds());

    label id = findIndex(idList, idName);

    if(id == -1)
    {
        FatalErrorIn("testInitialise::setInitialConfiguration()")
            << "Cannot find molecule id: " << idName << nl << "in idList."
            << exit(FatalError);
    }

    scalar massI = readScalar(propsDict_.lookup("mass"));

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
    
    scalar DX = Lx / Nx;
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

    
    // set exact number molecules?
    
    label N = 0;
    
    if (propsDict_.found("N"))
    {
        N = readLabel(propsDict_.lookup("N"));
        
        if(positions.size() < N)
        {
            FatalErrorIn("testInitialise::setInitialConfiguration()")
                << "Number of molecules in lattice  = " << positions.size()
                << ", number of molecules to allow are lower, N = " << N
                << exit(FatalError);            

        }
    }
    else
    {
        N = positions.size();
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
        
        vector v = equipartitionLinearVelocity(temperature, massI);
//         Info << "v = " << v << endl;
        
        v.z()=0.0;
        d.add(mag(v));
        
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
                frozen,
                v
            );
            
            nMolsInserted++;
        }
    }
    
    label nToDel = nMolsInserted - N;
    
    distributePoints randomBox
    (
        bb,
        cloud_.rndGen()
    );    
    
    // Delete excess molecules 
    
    if(nToDel > 0)
    {
        List<vector> molPositions(nToDel, vector::zero);
        
        forAll(molPositions, i)
        {
            molPositions[i] = randomBox.randomPoint();
        }
        
        DynamicList<agent*> molsToDel;
        DynamicList<label> chosenIds(0);
        
        forAll(molPositions, j)
        {   
            DynamicList<agent*> molsToDelTemp;
            
            const vector& rJ = molPositions[j];
            
            scalar deltaR = GREAT;        
            
            label tNI = 0;
            label chosenI = -1;
            
            IDLList<agent>::iterator mol(cloud_.begin());
            
            for
            (
                mol = cloud_.begin();
                mol != cloud_.end();
                ++mol
            )
            {
                if(mol().id() == id)
                {
                    scalar magRIJ = mag(rJ - mol().position());
                    agent* molI = &mol();
                    
                    if(magRIJ < deltaR)
                    {
                        if(findIndex(chosenIds, tNI) == -1)
                        {
                            deltaR = magRIJ;
                            
                            molsToDelTemp.clear();
                            
                            molsToDelTemp.append(molI);   
                            chosenI = tNI;
                        }
                    }
                }
                
                tNI++;
            }
            
            molsToDelTemp.shrink();
            
            if(chosenI != -1)
            {
                molsToDel.append(molsToDelTemp[0]);
                chosenIds.append(chosenI);
            }
        }  
        
        molsToDel.shrink();
        
        forAll(molsToDel, m)
        {
            cloud_.deleteParticle(*molsToDel[m]);
        }   
    }
    
    label finalSize = cloud_.size();

    nAgentsAdded_ = finalSize - initialSize;

    if (Pstream::parRun())
    {
        reduce(nAgentsAdded_, sumOp<label>());
    }

    Info << tab << " molecules added: " << nAgentsAdded_ << endl;
 
    
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

void testInitialise::setBoundBox
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
