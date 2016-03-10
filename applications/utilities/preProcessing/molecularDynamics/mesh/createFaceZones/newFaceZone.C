/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class
    newFaceZone

Description

\*----------------------------------------------------------------------------*/

#include "newFaceZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //




// Construct from IOdictionary and mesh
newFaceZone::newFaceZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    faces_()
{

    const word zoneName = dict.lookup("zoneName");
    regionName_ = zoneName;
        
    const faceZoneMesh& faceZones = mesh.faceZones();
    const label& regionId = faceZones.findZoneID(zoneName);
    
    if(regionId != -1)
    {
        FatalErrorIn("createFaceZones")
            << "FaceZone: " << zoneName << " exists on the mesh."
            << nl << " Check: "
            << "zoneDict"
            << nl << "Solve this by executing:"
            << nl << "> rm constant/polyMesh/*Zones"
            << nl << "> blockMesh"
            << nl << "> createFaceZones"            
            << exit(FatalError);
    }    
    
    
    if(dict.found("writeFaceSet"))
    {
        writeFaceSets_ = Switch(dict.lookup("writeFaceSet"));
    }

    const word option(dict.lookup("option"));
    option_ = option;
    
    setZone();    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

newFaceZone::~newFaceZone()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void newFaceZone::setZone()
{
    DynamicList<label> faces(0);    
    
    if(option_ == "boundBox")
    {

        vector startPoint = dict_.lookup("startPoint");
        vector endPoint = dict_.lookup("endPoint");
    
        boundedBox bb;
        bb.resetBoundedBox(startPoint, endPoint);
        
        for (label i = 0; i < mesh_.nFaces(); i++)
        {
            const vector& faceCentreI = mesh_.faceCentres()[i];        
            
            bool acceptedFace = false;
            
            if(bb.contains(faceCentreI))
            {
                faces.append(i);
                acceptedFace = true;
            }
            
            // further step of refinement can go here
            // Example if vertices of face + face centre amount to
            // > 50% inside boundbox, accept the face in the faceZone
            // FUTURE WORK
        }
    }
    
    if(option_ == "pointToPoint")
    {
        // not yet implemented 
        vector startPoint = dict_.lookup("startPoint");
        vector endPoint = dict_.lookup("endPoint");        
        vector rES = endPoint - startPoint;
       
        const vectorField& faceCentres = mesh_.faceCentres();

        forAll(faceCentres, f)
        {
            const vector& fC = faceCentres[f];

            scalar test1 = (fC - startPoint) & rES;
            scalar test2 = (endPoint - fC) & rES;

            if((test1 >= 0) && (test2 >= 0))
            {
                faces.append(f);
            }
        }
    }
        
    faces.shrink();

    faces_.setSize(faces.size());

    forAll(faces, i)
    {
        faces_[i] = faces[i];
    }     
    
//     Info << "-> number of faces: "<< faces_.size() << endl; 
    
    // Test can be placed here to check faces all contain the 
    // same normal to check if there are any bends on the plane
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
