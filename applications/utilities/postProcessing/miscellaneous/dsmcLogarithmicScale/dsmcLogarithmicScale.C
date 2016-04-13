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
   
    DynamicList<scalar> newTimeDyn(0);
    DynamicList<scalar> newPropertyDyn(0);
    
    label firstCounter = 0;
    
    //forAll(t, i)
    for(label r=0; r<times.size()-1; r++) 
    {
        scalar t1 = times[r];
        scalar t2 = times[r+1];
        
        label counter = 0;
        
        label x = label( (((t2-t1)/10)/deltaT) );

        forAll(time, i)
        {
            scalar t = time[i];
        
            if ( (t >= t1) && (t <= t2) )
            {
                if(firstCounter == 0)
                {
                    newTimeDyn.append(t);
                    newPropertyDyn.append(property[i]);
                    firstCounter++;
                }
                else
                {
                    counter++;
                    
                    if(counter == x)
                    {
                        newTimeDyn.append(t);
                        newPropertyDyn.append(property[i]);
                        counter = 0;
                    }                    
                }
            }
        }
    }
    
    fileName pathName(runTime.path());
    
    scalarField newTime;
    scalarField newProperty;

    newTime.transfer(newTimeDyn.shrink());
    newProperty.transfer(newPropertyDyn.shrink());
    
    writeTimeData
    (
        pathName,
        outputFileName,
        newTime,
        newProperty
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
