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
    To create a series of interfaces from any two touching zones.
    

Description


\*---------------------------------------------------------------------------*/

//#include "buildOverlap.H"

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"

#include "IFstream.H"

#include <iostream>
#include <iomanip>
#include <fstream>

#include "writeTimeData.H"

using namespace Foam;


// - Main program

int main(int argc, char *argv[])
{
//     argList::validArgs.append("patch");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

   Info << nl << "Reading dictionary" << endl;
   
   #   include "readDict.H"
   
    scalarField s(fftpoints, 0.0);
   
    scalarField taus(fftpoints,0);

    forAll(taus, i)
    {
        label t = i*tau;

        scalar summation = 0.0;

        label samples = 0;
        
        forAll(R, j)
        {
            if( (j+t) < ROWS )
            {
                summation += R[j]*R[j+t];
                samples++;
            }
        }
        
        summation /= samples;

        s[i] = summation;
        taus[i] = t;
    }
    
    fileName pathName(runTime.path());

    writeTimeData
    (
        pathName,
        outputFileName,
        taus,
        s
    );

    
    
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
