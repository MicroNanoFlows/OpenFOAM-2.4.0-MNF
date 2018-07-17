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
#include "molecules.H"
#include "allConfigurations.H"
#include "IFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "setTime.H"
#   include "createTime.H"
// #   include "createMesh.H"
    
//     Foam::fvMesh mesh
//     (
//         Foam::IOobject
//         (
//             Foam::fvMesh::defaultRegion,
//             "",
// //             runTime,
//             Foam::IOobject::NO_READ
//         )
//     );
    
    IOdictionary dict
    (
        IOobject
        (
            "in.spawn",
            "", //             mesh.time().system(),
            runTime, // mesh
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    
    molecules cloud
    (
//         mesh,
        dict
    );    

    allConfigurations confs
    (
//         mesh,
        cloud,
        dict
    );
    
    confs.spawn();

    cloud.write();

//     Info << nl << "ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl; 

    Info << nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
