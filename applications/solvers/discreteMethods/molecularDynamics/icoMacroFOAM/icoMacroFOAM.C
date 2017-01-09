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
    icoMacroFOAM

Description
    Incompressible macro solver for the Internal-flow Multiscale Method

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "icoMacroModel.H"

using namespace Foam;

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    reducedUnits redUnits(runTime, mesh);

    dictionary immDict =
    IOdictionary
    (
        IOobject
        (
            "immDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    icoMacroModel macro(runTime, redUnits, immDict);

    Info << "\nStarting time loop\n" << endl;
    
    
    for(label n = 0; n < macro.nIter(); n++)
    {
        // set and get forces 
        List<vector> forces = macro.setAndGetForces();
        
        // send forces to LAMMPS HERE
        
        //  RUN MD simulation HERE
//         macro.pseudoMD();
        List<scalar> mDot(macro.nMicro(),0.0);
        
        // get mass flow rates from MD
        // and set them in icoMacroModel 
        // change from LAMMPS units to SI
        
        macro.massFlowRates()=mDot;
        
        macro.solve();
        
        macro.write();
    }
    
    
    Info << "End\n" << endl;

    return 0;
}












