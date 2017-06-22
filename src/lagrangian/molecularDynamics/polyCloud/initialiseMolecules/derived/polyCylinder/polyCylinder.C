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

#include "polyCylinder.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyCylinder, 0);

addToRunTimeSelectionTable(polyConfiguration, polyCylinder, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyCylinder::polyCylinder
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
//     const word& name
)
:
    polyConfiguration(molCloud, dict/*, name*/)
//     propsDict_(dict.subDict(typeName + "Properties"))
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyCylinder::~polyCylinder()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polyCylinder::setInitialConfiguration()
{
    label initialSize = molCloud_.size();

    Info << nl << "Building quick simple lattice: " << endl;

    const scalar temperature(readScalar(mdInitialiseDict_.lookup("temperature")));

    const vector bulkVelocity(mdInitialiseDict_.lookup("bulkVelocity"));
    
    vector startPoint(mdInitialiseDict_.lookup("startPoint"));
    vector endPoint(mdInitialiseDict_.lookup("endPoint"));
    
    vector unitVector = (endPoint - startPoint)/mag(endPoint - startPoint);
    
    vector unitVectorX(mdInitialiseDict_.lookup("unitVectorX"));
    vector unitVectorY(mdInitialiseDict_.lookup("unitVectorY"));        
    
    
    scalar R(readScalar(mdInitialiseDict_.lookup("radius")));

    bool frozen = false;

    if (mdInitialiseDict_.found("frozen"))
    {
        frozen = Switch(mdInitialiseDict_.lookup("frozen"));
    }

    bool tethered = false;

    if (mdInitialiseDict_.found("tethered"))
    {
        tethered = Switch(mdInitialiseDict_.lookup("tethered"));
    }    
    
    const word molIdName(mdInitialiseDict_.lookup("molId")); 
    const List<word>& idList(molCloud_.cP().molIds());

    label molId = findIndex(idList, molIdName);

    if(molId == -1)
    {
        FatalErrorIn("polyCylinder::setInitialConfiguration()")
            << "Cannot find molecule id: " << molIdName << nl << "in idList."
            << exit(FatalError);
    }
    
    scalar massI = molCloud_.cP().mass(molId);
  
    
    boundedBox bb;
    
    
    vector rMin = startPoint - unitVectorX*R - unitVectorY*R;
    vector rMax = endPoint + unitVectorX*R + unitVectorY*R;
    
    bb.resetBoundedBox(rMin, rMax);

    scalar s(readScalar(mdInitialiseDict_.lookup("unitCellSize")));

    label nX = (bb.span().x()/s) + 1;
    label nY = (bb.span().y()/s) + 1;
    label nZ = (bb.span().z()/s) + 1;
    
    // basis points for bcc lattice which includes 2 sites
    label nAtoms= 0;

    DynamicList<vector> positions;
    
    for (label k = 0; k < nX; k++)
    {
        for (label j = 0; j < nY; j++)
        {
            for (label i = 0; i < nZ; i++)
            {
                vector p1 = vector(1, 0, 0)*k*s + vector(0, 1, 0)*j*s + vector(0, 0, 1)*i*s + bb.min();
                vector p2 = vector(1, 0, 0)*k*s + vector(1, 0, 0)*s*0.5 + 
                            vector(0, 1, 0)*j*s + vector(0, 1, 0)*s*0.5 +
                            vector(0, 0, 1)*i*s + vector(0, 0, 1)*s*0.5
                            + bb.min();
                            
                if(bb.contains(p1))
                {
                    positions.append(p1);
                    nAtoms++;
                }
                
                if(bb.contains(p2))
                {
                    positions.append(p2);
                    nAtoms++;
                }
            }
        }
    }
    
    positions.shrink();
    

    
    // reduce 
    
    DynamicList<vector> cutPositions;
    
    forAll(positions, i)
    {
        const vector& rI = positions[i];
        vector rIS = rI - startPoint;
        
        vector rD = rIS - (rIS & unitVector)*unitVector;
        
        if(mag(rD) <= R)
        {
            cutPositions.append(rI);
        }
    }
    
    positions.clear();
    
    positions.transfer(cutPositions);

    Info << nl << " No of sites found = " << positions.size() << endl;
    
    // set exact number molecules?
    
    label N = 0;
    

    N = readLabel(mdInitialiseDict_.lookup("N"));    
    
    if(positions.size() < N)
    {
        FatalErrorIn("polyBCC::setInitialConfiguration()")
            << "Number of molecules in lattice  = " << positions.size()
            << ", number of molecules to allow are lower, N = " << N
            << exit(FatalError);            

    }
    
    Info << "Target number of molecules to insert = " << N << endl;
    
    label nMolsInserted = 0;
    
    {
        Info << "Using randomness to select molecules" << endl;
        
        //label nToDel = positions.size()-N;
        
        // use probability to delete to get close to estimate
        
        scalar p = scalar(N)/scalar(positions.size());
        
        Info << "probability to accept molecules = " << p << endl;
        
        DynamicList<vector> rejectedPositions;
        
        
        // insert molecules
        forAll(positions, i)
        {
            if(molCloud_.rndGen().sample01<scalar>() <= p)
            {
                
                if(nMolsInserted < N)
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
                        insertMolecule
                        (
                            positions[i],
                            cell,
                            tetFace,
                            tetPt,
                            molId,
                            tethered,
                            frozen,
                            temperature,
                            bulkVelocity
                        );
                        
                        nMolsInserted++;
                    }
                }
            }
            else
            {
                rejectedPositions.append(positions[i]);
            }
        }
        
        Info << "Nmols inserted = " << nMolsInserted << endl;
        
        
        // add missing ones:
        forAll(rejectedPositions, i)
        {
            if(nMolsInserted < N)
            {
                label cell = -1;
                label tetFace = -1;
                label tetPt = -1;
                
                mesh_.findCellFacePt
                (
                    rejectedPositions[i],
                    cell,
                    tetFace,
                    tetPt
                );
                
                if(cell != -1)
                {
                    insertMolecule
                    (
                        rejectedPositions[i],
                        cell,
                        tetFace,
                        tetPt,
                        molId,
                        tethered,
                        frozen,
                        temperature,
                        bulkVelocity
                    );
                    
                    nMolsInserted++;
                }            
            }
        }
    }
    
    scalar V = bb.volume();

    scalar M = nMolsInserted*massI;
    
    scalar rhoM = M/V;
    
    Info << "Estimate of mass density, rhoM (RU) = " << rhoM 
        << ", SI = " << rhoM*molCloud_.redUnits().refMassDensity()
        << endl;


    label finalSize = molCloud_.size();

    nMolsAdded_ = finalSize - initialSize;

    if (Pstream::parRun())
    {
        reduce(nMolsAdded_, sumOp<label>());
    }

    Info << tab << " molecules added: " << nMolsAdded_ << endl;
    
}



} // End namespace Foam

// ************************************************************************* //
