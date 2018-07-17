/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "unitCell.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(unitCell, 0);

addToRunTimeSelectionTable(configuration, unitCell, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
unitCell::unitCell
(
    molecules& cloud,
    const dictionary& dict
)
:
    configuration(cloud, dict)
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

unitCell::~unitCell()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void unitCell::spawn()
{
    word fileName(dict_.lookup("fileName"));
    
    IFstream data(fileName);
    
    label nAtoms;
    
    scalar x = 0.0;
    scalar y = 0.0;
    scalar z = 0.0;
    scalar q = 0.0;
    word type;    
     
    Info << nl << "reading unit cell from filename " << fileName << endl;
    
    // read in unit cell structure consisting of: type x y z q
    
    DynamicList<vector> unitCellPositions;
    DynamicList<label> unitCellTypes;
    DynamicList<scalar> unitCellCharges;
    
    data >> nAtoms;
    
    for(label i=0; i<nAtoms; i++)
    {
        data >> type >> x >> y >> z >> q;
        
        Info<< type  << " " 
            << x << " " 
            << y << " " 
            << z << " " 
            << q << " " 
            << endl;
            
        unitCellPositions.append(vector(x, y, z));
        unitCellCharges.append(q);
        
       
        if(findIndex(cloud_.typeNames(), type) == -1)
        {
            FatalErrorIn("unitCell::spawn()")
                << "Cannot find atom type: " << type
                << nl << "in types list."
                << exit(FatalError);            
        }
        
        label typeId = cloud_.getType(type);
        
        unitCellTypes.append(typeId);
    }
    
    label nX(readScalar(dict_.lookup("nX")));
    label nY(readScalar(dict_.lookup("nY")));
    label nZ(readScalar(dict_.lookup("nZ")));
    
    scalar Lx(readScalar(dict_.lookup("Lx")));
    scalar Ly(readScalar(dict_.lookup("Ly")));
    scalar Lz(readScalar(dict_.lookup("Lz")));
    
    label nSites = nAtoms;
    

    // Bounds box 

    boundsBox bb;
    
    if(dict_.found("boundsBox"))
    {
        setBoundsBox(dict_, bb, "boundsBox");
    }
    else
    {
        bb = cloud_.bMesh();
    }
    
    // some checks 
    
    scalar LxBox = bb.span().x();
    scalar LyBox = bb.span().y();
    scalar LzBox = bb.span().z();
    
    if(Lx*nX > LxBox)
    {
        FatalErrorIn("unitCell::spawn()")
            << "Your box is too small for nX = "<<nX<<" replications. Try:"
            << nl
            << "a) nX = " << floor(LxBox/Lx) << " or, "
            << "b) making box bigger = " << LxBox 
            << exit(FatalError);       
    }
    
    if(Ly*nY > LyBox)
    {
        FatalErrorIn("unitCell::spawn()")
            << "Your box is too small for nY = "<<nY<<" replications. Try:"
            << nl
            << "a) nY = " << floor(LyBox/Ly) << " or, "
            << "b) making box bigger = " << LyBox 
            << exit(FatalError);      
    }
    
    if(Lz*nZ > LzBox)
    {
        FatalErrorIn("unitCell::spawn()")
            << "Your box is too small for nZ = "<<nZ<<" replications. Try:"
            << nl
            << "a) nZ = " << floor(LzBox/Lz) << " or, "
            << "b) making box bigger = " << LzBox 
            << exit(FatalError);  
    }
    

    label noOfAtoms= 0;
    
    DynamicList<vector> positions;
    DynamicList<label> types;
    DynamicList<scalar> charges;
    
    for (label k = 0; k < nX; k++)
    {
        for (label j = 0; j < nY; j++)
        {
            for (label i = 0; i < nZ; i++)
            {
                for (label iS = 0; iS < nSites; iS++)
                {
                    vector rI = unitCellPositions[iS];
                    
                    scalar xS = rI.x() + Lx*k + bb.min().x();
                    scalar yS = rI.y() + Ly*j + bb.min().y();
                    scalar zS = rI.z() + Lz*i + bb.min().z();
                                        
                    vector pos = vector(xS, yS, zS);

                    
                    if(bb.contains(pos))
                    {
                        positions.append(pos);
                        types.append(unitCellTypes[iS]);
                        charges.append(unitCellCharges[iS]);
                        noOfAtoms++;
                    }
                    else
                    {
                        // error
                    }
                }
            }
        }
    }
    

    Info << nl << " Number of sites found = " << positions.size() << endl;

    // check for overlaps here (to add)
    
    

    forAll(positions, i)
    {
        cloud_.positions().append(positions[i]);
        cloud_.types().append(types[i]);
        cloud_.charges().append(charges[i]);
    }
    
}

void unitCell::setBoundsBox
(
    const dictionary& propsDict,
    boundsBox& bb,
    const word& name 
)
{
    const dictionary& dict(propsDict.subDict(name));
    
    vector startPoint = dict.lookup("startPoint");
    vector endPoint = dict.lookup("endPoint");

    bb.resetBoundedBox(startPoint, endPoint);
}




} // End namespace Foam

// ************************************************************************* //
