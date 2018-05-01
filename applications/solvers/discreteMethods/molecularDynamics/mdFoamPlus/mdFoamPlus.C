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

Application
    mdPolyFoam

Description
    Molecular dynamics solver for fluid dynamics -- polyatomics only

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#ifdef USE_MUI // included if the switch -DUSE_MUI included during compilation.
  #include "fvCoupling.H"
#endif
#include "mdPoly.H"
#ifdef USE_MUI // included if the switch -DUSE_MUI included during compilation.
    #include "mui.h"
#endif

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createRandom.H"

#ifdef USE_MUI // included if the switch -DUSE_MUI included during compilation.
  #   include "createCouplingData.H"

  if (args.cplRunControl().cplRun())
  {
    #   include "createCouplings.H"
  }
#endif
    
    reducedUnits rU(runTime, mesh);

    constantMoleculeProperties cP (mesh, rU);

    polyMoleculeCloud *molecules;
        
#ifdef USE_MUI
    if (args.cplRunControl().cplRun())
    {
        molecules = new polyMoleculeCloud
        (
            runTime,
            mesh,
            rU,
            cP,
            rndGen,
            oneDInterfaces,
            twoDInterfaces,
            threeDInterfaces
        );
    }
    else
    {
        molecules = new polyMoleculeCloud
        (
            runTime,
            mesh,
            rU,
            cP,
            rndGen
        );
    }
#else
    molecules = new polyMoleculeCloud
    (
        runTime,
        mesh,
        rU,
        cP,
        rndGen
    );
#endif

    Info << "\nStarting time loop\n" << endl;

//     clockTimer evolveTimer(runTime, "evolve", true);

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << endl;

        molecules->clock().startClock();

        molecules->evolve();

        molecules->clock().stopClock();

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info << "End\n" << endl;

#ifdef USE_MUI // included if the switch -DUSE_MUI included during compilation.
    #include "deleteCouplings.H"
#endif

    return 0;
}
