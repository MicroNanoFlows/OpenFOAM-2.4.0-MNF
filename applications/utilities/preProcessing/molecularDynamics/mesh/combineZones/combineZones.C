/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 1991-2005 OpenCFD Ltd.
    \\/      M anipulation   |
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

Application
    To combine a series of zones together into one zone.
    AND
    To subtract the combined zone from the global zone.
    

Description


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"

#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"

using namespace Foam;

//-sets
#include "createCellSet.H"
//-zones
#include "createCellZone.H"
#include "createFaceZone.H"
//-conversionsToFaceZones
#include "cellsToFaceZone.H"

// - Main program
int main(int argc, char *argv[])
{

#	include "setRootCase.H"
#	include "createTime.H"
#	include "createPolyMesh.H"

    Info << nl << "Reading combineZonesDict" << nl << endl;

#	include "readCombineZonesDict.H"

    //- global region

    List<label> globalRegion(mesh.nCells());

    for(label i = 0; i < globalRegion.size(); i++)
    {
        globalRegion[i] = i;
    }

    forAll(regionsToCombineList, rC)
    {
        DynamicList<label> addedRegions(0);

        Info << rC<< " combine zone: " << addedRegionsNames[rC] << endl;

        forAll(regionZoneNames[rC], rZ)
        {
            const word& cellZoneName = regionZoneNames[rC][rZ];
            const label& cellZoneId = cellZones.findZoneID(cellZoneName);
            const labelList& cellList = cellZones[cellZoneId];
    
            forAll(cellList, c)
            {
                const label& cellI = cellList[c];
    
                if(findIndex(addedRegions, cellI) == -1)
                {
                    addedRegions.append(cellI);
                }
            }
        }
    
        addedRegions.shrink();
    
        createCellZone(mesh, addedRegions, addedRegionsNames[rC]);
    
        // create sets
        if(writeCellSets[rC])
        {
            createCellSet(mesh, addedRegions, addedRegionsNames[rC]);
        }
        
        // visualisation of region
        cellsToFaceZone(mesh, addedRegions, addedRegionsNames[rC]);
    
    
        // - subtract added region from the global domain
    
        if(writeSubtractedRegions[rC])
        {
            DynamicList<label> subtractedRegion(0);
        
            forAll(globalRegion, c)
            {
                const label& cellI = globalRegion[c];
        
                if(findIndex(addedRegions, cellI) == -1)
                {
                    subtractedRegion.append(cellI);
                }
            }
        
            subtractedRegion.shrink();
            
            createCellZone(mesh, subtractedRegion, subtractedRegionsNames[rC]);
    
            // create sets
            if(writeCellSets[rC])
            {
                createCellSet(mesh, subtractedRegion, subtractedRegionsNames[rC]);
            }
        
            // visualisation of region
            cellsToFaceZone(mesh, subtractedRegion, subtractedRegionsNames[rC]);
        }
    }

    if (!mesh.write())     // - requrired for zone writing
    {
            FatalErrorIn(args.executable())
                << "Failed writing cellZones."
                << exit(FatalError);
    }

    Info << nl << "ClockTime = " << runTime.elapsedClockTime() 
               << " s" << nl << endl; 

    Info << nl << "end\n" << endl;

    return 0;
}

// ************************************************************************* //
