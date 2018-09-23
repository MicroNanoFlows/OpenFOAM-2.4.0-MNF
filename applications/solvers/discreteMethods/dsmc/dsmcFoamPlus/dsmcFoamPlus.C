/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    dsmcFoam

Description
    Direct simulation Monte Carlo (DSMC) solver for 3D, transient, multi-
    species flows

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
//     #include "createDynamicFvMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Constructing dsmcCloud " << endl;

    dsmcCloud dsmc(runTime, "dsmc", mesh); 

    Info<< "\nStarting time loop\n" << endl;
    
    label infoCounter = 0;

    while (runTime.loop())
    {          
        infoCounter++;
        
        if(infoCounter >= dsmc.nTerminalOutputs())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
        }

        dsmc.evolve();

        if(infoCounter >= dsmc.nTerminalOutputs())
        {
            dsmc.info();   
        }

        runTime.write();

        if(infoCounter >= dsmc.nTerminalOutputs())
        {
            Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
                
            infoCounter = 0;
        }
        
        dsmc.loadBalanceCheck();
            
//         scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();
// 
//         mesh.update();
// 
//         if (mesh.changing())
//         {
//             Info<< "Execution time for mesh.update() = "
//                 << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
//                 << " s" << endl;
//         }
    }

    Info<< "End\n" << endl;
    
    dsmc.loadBalance();

    return(0);
}


// ************************************************************************* //
