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
    Deletes molecules based on an input model.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mdPoly.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"

    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"
#   include "createRandom.H"
    

    reducedUnits rU(runTime, mesh);

    constantMoleculeProperties cP (mesh, rU);
        
    polyMoleculeCloud molecules
    (
        runTime,
        mesh,
        rU,
        cP,
        rndGen,
        "delete",
        false
    );    
    
    IOstream::defaultPrecision(15);

    runTime++;

    Info << nl << "Writing fields." << endl;

    if (!mesh.write())
    {
        FatalErrorIn(args.executable())
            << "Failed writing moleculeCloud."
            << exit(FatalError);
    }
    
    Info << nl << "ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl; 

    Info << nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
