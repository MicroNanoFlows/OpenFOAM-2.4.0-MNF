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

#include "polyBCC.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyBCC, 0);

addToRunTimeSelectionTable(polyConfiguration, polyBCC, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyBCC::polyBCC
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyConfiguration(molCloud, dict)
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyBCC::~polyBCC()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polyBCC::setInitialConfiguration()
{
    label initialSize = molCloud_.size();

    Info << nl << "BCC lattice " << endl;

    bool SIUnits = false;

    if (mdInitialiseDict_.found("SIUnits"))
    {
        SIUnits = Switch(mdInitialiseDict_.lookup("SIUnits"));
    }

    scalar temperature(readScalar(mdInitialiseDict_.lookup("temperature")));

    if(SIUnits)
    {
        if(temperature != 0)
        {
            temperature /= molCloud_.redUnits().refTemp();
        }
    }

    vector bulkVelocity(mdInitialiseDict_.lookup("bulkVelocity"));

    if(SIUnits)
    {
        if(bulkVelocity[0] != 0)
        {
            bulkVelocity[0] /= molCloud_.redUnits().refVelocity();
        }

        if(bulkVelocity[1] != 0)
        {
            bulkVelocity[1] /= molCloud_.redUnits().refVelocity();
        }

        if(bulkVelocity[2] != 0)
        {
            bulkVelocity[2] /= molCloud_.redUnits().refVelocity();
        }
    }

    const word molIdName(mdInitialiseDict_.lookup("molId")); 
    const List<word>& idList(molCloud_.cP().molIds());

    label molId = findIndex(idList, molIdName);

    if(molId == -1)
    {
        FatalErrorIn("polyBCC::setInitialConfiguration()")
            << "Cannot find molecule id: " << molIdName << nl << "in idList."
            << exit(FatalError);
    }

    scalar massI = molCloud_.cP().mass(molId);

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
    
    bool multispecies = false;
    
    if(mdInitialiseDict_.found("multispecies"))
    {
        multispecies = Switch(mdInitialiseDict_.lookup("multispecies"));
    }
    
    // bounding box 

    boundedBox bb;
    
    setBoundBox(mdInitialiseDict_, bb, "boundBox");

    if(SIUnits)
    {
        bb.min()[0] /= molCloud_.redUnits().refLength();
        bb.min()[1] /= molCloud_.redUnits().refLength();
        bb.min()[2] /= molCloud_.redUnits().refLength();

        bb.max()[0] /= molCloud_.redUnits().refLength();
        bb.max()[1] /= molCloud_.redUnits().refLength();
        bb.max()[2] /= molCloud_.redUnits().refLength();
    }

    scalar s(readScalar(mdInitialiseDict_.lookup("unitCellSize")));

    label nX = static_cast<label>(ceil(bb.span().x()/s)) + 1;
    label nY = static_cast<label>(ceil(bb.span().y()/s)) + 1;
    label nZ = static_cast<label>(ceil(bb.span().z()/s)) + 1;
    
    // basis points for bcc lattice which includes 2 sites
    label nAtoms= 0;

    DynamicList<vector> positions;
    
    for (label k_i = 0; k_i < nX; k_i++)
    {
        for (label j_i = 0; j_i < nY; j_i++)
        {
            for (label i_i = 0; i_i < nZ; i_i++)
            {
                scalar i = static_cast<scalar>(i_i);
                scalar j = static_cast<scalar>(j_i);
                scalar k = static_cast<scalar>(k_i);

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
    

    Info << nl << " No of sites found = " << positions.size() << endl;

    // set exact number molecules?
    
    label N = 0;
    
    if (mdInitialiseDict_.found("N"))
    {
        N = readLabel(mdInitialiseDict_.lookup("N"));
        
        if(positions.size() < N)
        {
            FatalErrorIn("polyBCC::setInitialConfiguration()")
                << "Number of molecules in lattice  = " << positions.size()
                << ", number of molecules to allow are lower, N = " << N
                << exit(FatalError);            

        }
        
        Info << "Target number of molecules to insert = " << N << endl;
    }
    else
    {
        N = positions.size();
    }
    
    label nMolsInserted = 0;
    
    if(N == positions.size())
    {
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
        Info << "Using randomness to select molecules" << endl;
        
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
    
    
    //- new - delete existing molecules and replace with other species
    if(multispecies)
    {
        distributePoints randomBox
        (
          bb,
          molCloud_.rndGen()
        );
        
        selectIds ids
        (
           molCloud_.cP(),
           mdInitialiseDict_
        );
        
        List<label> molIds = ids.molIds();
        List<label> soluteN = List<label>(mdInitialiseDict_.lookup("soluteN"));   
        
        forAll(molIds, i)
        {
            List<vector> molPositions(soluteN[i], vector::zero);
            
            forAll(molPositions, j)
            {
                molPositions[j] = randomBox.randomPoint();
            }
            
            DynamicList<polyMolecule*> molsToDel;
            DynamicList<label> chosenIds(0);
            
            forAll(molPositions, j)
            {   
                DynamicList<polyMolecule*> molsToDelTemp;
                
                const vector& rJ = molPositions[j];
                
                scalar deltaR = GREAT;        
                
                label tNI = 0;
                label chosenI = -1;
                
                IDLList<polyMolecule>::iterator mol(molCloud_.begin());
                
                for
                (
                     mol = molCloud_.begin();
                     mol != molCloud_.end();
                     ++mol
                )
                {
                    if(mol().id() == molId)
                    {
                        scalar magRIJ = mag(rJ - mol().position());
                        polyMolecule* molI = &mol();
                        
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
            
            DynamicList<vector> newPositions;
            
            forAll(molsToDel, m)
            {
                newPositions.append(molsToDel[m]->position());
                molCloud_.deleteParticle(*molsToDel[m]);
            }            
            
            newPositions.shrink();
            label nExchanged = 0;
            
            forAll(newPositions, j)
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
                        molIds[i],
                        tethered,
                        frozen,
                        temperature,
                        bulkVelocity
                    );
                    
                    nExchanged++;
                }
            }
            
            Info<< "Mol id = " << molIds[i] 
                << " Number of molecules exchanged = " << nExchanged
                << endl;
        }
    }
}

void polyBCC::setBoundBox
(
    const dictionary& propsDict,
    boundedBox& bb,
    const word& name 
)
{
    const dictionary& dict(propsDict.subDict(name));
    
    vector startPoint = dict.lookup("startPoint");
    vector endPoint = dict.lookup("endPoint");

    bool SIUnits = false;

    if (dict.found("SIUnits"))
    {
        SIUnits = Switch(mdInitialiseDict_.lookup("SIUnits"));
    }

    if(SIUnits)
    {
        startPoint /= molCloud_.redUnits().refLength();
        endPoint /= molCloud_.redUnits().refLength();
    }

    bb.resetBoundedBox(startPoint, endPoint);
}




} // End namespace Foam

// ************************************************************************* //
