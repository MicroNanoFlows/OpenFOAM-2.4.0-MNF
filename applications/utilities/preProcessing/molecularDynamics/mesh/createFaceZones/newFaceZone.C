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

void newFaceZone::checkBoundBox
(
    boundBox& b,
    const vector& startPoint,
    const vector& endPoint
)
{
    vector& vMin = b.min();
    vector& vMax = b.max();

    if(startPoint.x() < endPoint.x())
    {
        vMin.x() = startPoint.x();
        vMax.x() = endPoint.x();
    }
    else
    {
        vMin.x() = endPoint.x();
        vMax.x() = startPoint.x();
    }
    if(startPoint.y() < endPoint.y())
    {
        vMin.y() = startPoint.y();
        vMax.y() = endPoint.y();
    }
    else
    {
        vMin.y() = endPoint.y();
        vMax.y() = startPoint.y();
    }
    if(startPoint.z() < endPoint.z())
    {
        vMin.z() = startPoint.z();
        vMax.z() = endPoint.z();
    }
    else
    {
        vMin.z() = endPoint.z();
        vMax.z() = startPoint.z();
    }
}       
    
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //




// Construct from IOdictionary and mesh
newFaceZone::newFaceZone
(
    const polyMesh& mesh,
    const word& zoneName,
    const vector& startPoint,
    const vector& endPoint
)
:
    mesh_(mesh),
    regionName_(zoneName),
    startPoint_(startPoint),
    endPoint_(endPoint),
    faces_()
{
    checkBoundBox
    (
        bb_,
        startPoint_,
        endPoint_
    );
    
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
    
    for (label i = 0; i < mesh_.nFaces(); i++)
    {
        const vector& faceCentreI = mesh_.faceCentres()[i];        
        
        bool acceptedFace = false;
        
        if(bb_.contains(faceCentreI))
        {
            faces.append(i);
            acceptedFace = true;
        }
        
        // further step of refinement can go here
        // Example if vertices of face + face centre amount to
        // > 50% inside boundbox, accept the face in the faceZone
        // FUTURE WORK
    }
    
    faces.shrink();

    faces_.setSize(faces.size());

    forAll(faces, i)
    {
        faces_[i] = faces[i];
    }    
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
