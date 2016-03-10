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
    mapTwoLagrangianFields

Description
    Copy an existing cloud of molecules, which has been obtained 
    from a different mesh/case in the 0 time-directory of the new case.

    In the next time-step directory (e.g. 0.005) have an existing cloud of 
    molecules which has already been mapped (e.g. using mapLagrangianFields 
    utility) to the current mesh/case.

    The utlity then checks that the molecule-to-cell addressing is correct only 
    on the new field (read in the 0 time directory), and combines the two
    Lagrangian fields together in the next time directory (e.g. 0.01).

    Not parallelised. Doesn't necessarily need to be, since it is not 
    computationally demanding.


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
        "NULL",
        false
    );    

    runTime++;

    label initialSize = molecules.size();

    {
        
        polyMoleculeCloud oldMolecules
        (
            runTime,
            mesh,
            rU,
            cP,
            rndGen,
            "mapping",
            false
        );    
        
        IDLList<polyMolecule>::iterator mol(oldMolecules.begin());

        for
        (
            mol = oldMolecules.begin();
            mol != oldMolecules.end();
            ++mol
        )
        {
            molecules.addParticle
            (
                new polyMolecule
                (
                    mol()
                )
            );
        }
    }
   
    runTime++;

    Info << nl << "Original no. of molecules: " << initialSize 
         << ", combined no. of molecules: " << molecules.size() 
         << endl;


    Info << nl << "Writing fields." << endl;

    IOstream::defaultPrecision(12);

    if (!mesh.write())
    {
        FatalErrorIn(args.executable())
            << "Failed writing atomisticMoleculeCloud."
            << exit(FatalError);
    }
    
    Info << nl << "ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl; 

    Info << nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
