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
    newCellZone

Description

\*----------------------------------------------------------------------------*/

#include "newCellZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void newCellZone::checkBoundBox
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
newCellZone::newCellZone
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
    cells_()
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

newCellZone::~newCellZone()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void newCellZone::setZone()
{
    DynamicList<label> cells(0);
    
    for (label i = 0; i < mesh_.nCells(); i++)
    {
        const vector& cellCentreI = mesh_.cellCentres()[i];        
        
        bool acceptedCell = false;
        
        if(bb_.contains(cellCentreI))
        {
            cells.append(i);
            acceptedCell = true;
        }
        
        // further step of refinement can go here
        
    }
    
    cells.shrink();

    cells_.setSize(cells.size());

    forAll(cells, i)
    {
        cells_[i] = cells[i];
    }    

    Info << "Target bound box: " 
         << bb_
         << endl;

    //- computing the size of the square

    scalar minX = GREAT;
    scalar maxX = 0.0;
    scalar minY = GREAT;
    scalar maxY = 0.0;
    scalar minZ = GREAT;
    scalar maxZ = 0.0;


    forAll(cells_, c)
    {
        const label& cellI = cells_[c];
        const labelList& pointList = mesh_.cellPoints()[cellI];

        forAll(pointList, p)
        {
            const label& pointI = pointList[p];
            const vector& vPointI = mesh_.points()[pointI];

            if(vPointI.x() < minX)
            {
                minX = vPointI.x();
            }

            if(vPointI.y() < minY)
            {
                minY = vPointI.y();
            }

            if(vPointI.z() < minZ)
            {
                minZ = vPointI.z();
            }

            if(vPointI.x() > maxX)
            {
                maxX = vPointI.x();
            }

            if(vPointI.y() > maxY)
            {
                maxY = vPointI.y();
            }

            if(vPointI.z() > maxZ)
            {
                maxZ = vPointI.z();
            }
        }

    }

    Info << "Dimensions of square: " 
         << mag(maxX - minX) << " x " 
         << mag(maxY - minY) << " x " 
         << mag(maxZ - minZ)
         << endl;
         
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
