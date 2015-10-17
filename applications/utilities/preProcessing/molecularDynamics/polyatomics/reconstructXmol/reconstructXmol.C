/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2008-2009 OpenCFD Ltd.
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

Description
    Reconstructs the .xmol files required for visualisation in VMD, paraView etc.

    For each time-directory the utlitity is used, the cloud is read in from the hard
    disc and the .xmol file uses the read in positions to write it out.

    Use examples:
    (for all time directories)
        > reconstrucXmol 

    (for one time directory, e.g. 12)
        > reconstrucXmol -time 12

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mdPoly.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

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
        "NULL",
        false
    );

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time: " << runTime.timeName() << endl;

        molecules.readNewField();

        fileName path(runTime.path()/runTime.timeName() +  + "/lagrangian" + "/polyMoleculeCloud.xmol");

        molecules.writeXYZ(path);

        runTime.write();
    }

    Info<< nl << "ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info << nl << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
