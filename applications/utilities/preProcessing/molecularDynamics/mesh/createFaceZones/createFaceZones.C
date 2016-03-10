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
    

Description


\*---------------------------------------------------------------------------*/

//#include "buildOverlap.H"

#include "argList.H"
#include "fvCFD.H"
#include "Time.H"
#include "polyMesh.H"

#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "newFaceZone.H"


using namespace Foam;
// #   include "backUpFunctions.H"
//-sets
#   include "createCellSet.H"
#   include "createFaceSet.H"
//-zones
#   include "createCellZone.H"
#   include "createFaceZone.H"

//-conversionsToFaceZones
#   include "cellsToFaceZone.H"

// - Main program

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"


    Info << nl << "Reading faceZonesDict" << nl<< endl;

#   include "readFaceZonesDict.H"

    forAll(surfacesList, sL)
    {
        Info << sL << ": Creating face zone: " << surfaceNames[sL] << endl;

        newFaceZone surfaceI
        (
            mesh,
            dictionaries[sL]
        );

        Info << " number of faces: " << surfaceI.faces().size() << nl<< endl;

        createFaceZone(mesh, surfaceI.faces(), surfaceI.name());
    
        // create sets
        if(surfaceI.writeFaceSets())
        {
            createFaceSet(mesh, surfaceI.faces(), surfaceI.name());
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
