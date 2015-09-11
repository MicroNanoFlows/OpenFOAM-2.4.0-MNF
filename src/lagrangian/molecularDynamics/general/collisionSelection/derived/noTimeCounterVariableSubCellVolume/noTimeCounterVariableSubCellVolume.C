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
    noTimeCounterVariableSubCellVolume

Description

\*----------------------------------------------------------------------------*/

#include "noTimeCounterVariableSubCellVolume.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noTimeCounterVariableSubCellVolume, 0);

addToRunTimeSelectionTable(collisionSelection, noTimeCounterVariableSubCellVolume, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

pointField noTimeCounterVariableSubCellVolume::cellPoints(const label& cell)
{
    const labelList& points = mesh_.cellPoints()[cell];

    pointField vectorPoints(points.size(), vector::zero);

    forAll(points, p)
    {
        vectorPoints[p] = mesh_.points()[points[p]];
    }

    return vectorPoints;
}

void noTimeCounterVariableSubCellVolume::checkSubCelling()
{


    // need to output these property in geometricFields to visualise level 
    // of skewness, errors etc in paraView

    // Test 1: check for ratio of cell-bounding box volume to actual cell volume.

    // Test 2: 

    for (label c = 0; c < mesh_.nCells(); c++)
    {
        const scalar& cellV = mesh_.cellVolumes()[c];

        boundBox bb(cellPoints(c), false);
        scalar bbVol = bb.span().x() * bb.span().y() * bb.span().z();

        scalar rV = cellV/bbVol;

        scalar nTot = bbVol/subCellVolume_;
        scalar rN = nTot/nSubCells_[c];

        scalar vSC = binWidths_[c].x()*binWidths_[c].y()*binWidths_[c].z();
        scalar rSC = subCellVolume_/vSC;

        Info << "cell " << c << nl
            << "Test 1 - cell volume ratios. " << nl 
            <<"cell volume = " << cellV
            << ", bb volume = " << bbVol
            << ", volume ratio = " << rV
            << nl
            << "Test 2 - no of sub cells. " << nl
            << "no. of required sub cells: " << nTot
            << ", no. of placed sub cells = " << nSubCells_[c]
            << ", ratio of error: " << rN
            << nl
            << "Test 3 - sub cell volume comparisons. " << nl
            << "target sub-cell volume: " << subCellVolume_
            << ", actual sub-cell volume: " << vSC
            << ", ratio of error: " << rSC
            << endl;
    }

    // check for number of sub-cells per cell on the mesh

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
noTimeCounterVariableSubCellVolume::noTimeCounterVariableSubCellVolume
(
    const polyMesh& mesh,
//     dsmcCloud& cloud
    const dictionary& dict
)
:
    collisionSelection(mesh/*, cloud*/, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    subCellVolume_(readScalar(propsDict_.lookup("subCellVolume"))),
    startPoints_(mesh.nCells(), vector::zero),
    nSlices_(mesh.nCells()),
    nSubCells_(mesh.nCells(), -1),
    binWidths_(mesh.nCells(), vector::zero)
{
    forAll(nSlices_, c)
    {
        nSlices_[c].setSize(3, 0);
    }

    if(subCellVolume_ <= 0)
    {
        Info << "error: choose positive sub-cell volume" << endl;
    }
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noTimeCounterVariableSubCellVolume::~noTimeCounterVariableSubCellVolume()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noTimeCounterVariableSubCellVolume::initialConfiguration()
{
//--TEMP
//     {
//         pointField p(8);
// 
//         p[0] = point(0, 0, 0);
//         p[1] = point(1, 1, 0);
//         p[2] = point(0, 2, 0);
//         p[3] = point(-1, 1, 0);
//         p[4] = point(0, 0, 2);
//         p[5] = point(1, 1, 2);
//         p[6] = point(0, 2, 2);
//         p[7] = point(-1, 1, 2);
// 
//         Info << "p" << p << endl;
//         boundBox bb(p, false);
// 
//         Info<< "bb min : " << bb.min() 
//             << ", bb max : " << bb.max() 
//             << ", span: " << bb.span()
//             << ", mid-point: " << bb.midpoint()
//             << endl;
//     }
//--
    // 1. Create boundbox around each cell and find the starting point
    //    i.e. the minimum point at the base of the bound box

    // 2. Identify the number of slices in the x,y z coordinates and their binWidths

    for (label c = 0; c < mesh_.nCells(); c++)
    {
        boundBox bb(cellPoints(c), false);

        startPoints_[c] = bb.min();

        scalar bbVol = bb.span().x() * bb.span().y() * bb.span().z();

        // total no of sub cells
        scalar nTot = bbVol/subCellVolume_;

        Info << "bound box, cC: " << bb.midpoint() 
            << ", start point: " << startPoints_[c]
            << ", bbVol: " << bbVol
            << ", total no of sub-cells = " << nTot 
            << endl;

        // test for no subcelling
        if(nTot <= 1.0)
        {
            nSlices_[c][0] = 1;
            nSlices_[c][1] = 1;
            nSlices_[c][2] = 1;
        }
        else
        {
            scalar deltaLCubic = Foam::pow(bbVol, (1.0/3.0) );
            scalar spacingCubic = deltaLCubic/Foam::pow(nTot, (1.0/3.0) );
    
            scalar nX = bb.span().x()/spacingCubic;
            scalar nY = bb.span().y()/spacingCubic;
            scalar nZ = bb.span().z()/spacingCubic;
    
            nSlices_[c][0] = label(nX+0.5);
            nSlices_[c][1] = label(nY+0.5);
            nSlices_[c][2] = label(nZ+0.5);
    
            // test for zero slices (1 is the bare minimum) and modfiy
            if(nSlices_[c][0] == 0)
            {
                nSlices_[c][0] = 1;
            }
            if(nSlices_[c][1] == 0)
            {
                nSlices_[c][1] = 1;
            }
            if(nSlices_[c][2] == 0)
            {
                nSlices_[c][2] = 1;
            }

            Info << "slices: nX = " << nSlices_[c][0]
                << ", nY = " << nSlices_[c][1]
                << ", nZ = " << nSlices_[c][2]
                << endl;
        }

        binWidths_[c].x() = bb.span().x()/scalar(nSlices_[c][0]);
        binWidths_[c].y() = bb.span().y()/scalar(nSlices_[c][1]);
        binWidths_[c].z() = bb.span().z()/scalar(nSlices_[c][2]);

        nSubCells_[c] = nSlices_[c][0]*nSlices_[c][1]*nSlices_[c][2];

        Info << "widths: nX = " << binWidths_[c].x()
            << ", nY = " << binWidths_[c].y()
            << ", nZ = " << binWidths_[c].z()
            << endl;
    }

    checkSubCelling();
}

void noTimeCounterVariableSubCellVolume::collide()
{
    //test for cell 0
    label c = 0;

    List<DynamicList<label> > subCells(nSubCells_[c]);

    List<vector> cellParcels(27);

    label counter = 0;

    for (label i = 0; i < nSlices_[c][0]; i++)
    {
        scalar x = i*binWidths_[c].x() + 0.5*binWidths_[c].x();

        for (label j = 0; j < nSlices_[c][1]; j++)
        {
            scalar y = j*binWidths_[c].y() + 0.5*binWidths_[c].y();

            for (label k = 0; k < nSlices_[c][2]; k++)
            {
                scalar z = k*binWidths_[c].z() + 0.5*binWidths_[c].z();
    
                cellParcels[counter] = startPoints_[c] + x*vector(1, 0, 0) 
                                        + y*vector(0, 1, 0) + z*vector(0, 0, 1);
                counter++;
            }
        }
    }

    Info << " nSubCells: " << nSubCells_[c] 
         << ", cell parcels: " << cellParcels
         << endl; 


    List<label> whichSubCell(cellParcels.size());

    forAll(cellParcels, i)
    {
        vector pS = cellParcels[i] - startPoints_[c];

        label nX = label((pS & vector(1, 0, 0))/binWidths_[c].x());

        label nY = label((pS & vector(0, 1, 0))/binWidths_[c].y());

        label nZ = label((pS & vector(0, 0, 1))/binWidths_[c].z());

        label subCell = nX + nY*nSlices_[c][0] + nZ*nSlices_[c][0]*nSlices_[c][1];

        Info << "position: "<< cellParcels[i]
            << ", nX = " << nX
            << ", nY = " << nY
            << ", nZ = " << nZ
            << ", subCell index = " << subCell
            << endl;

        subCells[subCell].append(i);

        whichSubCell[i] = subCell;
    }

    Info << "subCells: " << subCells << endl;
    Info << "whichSubCell: " << whichSubCell << endl;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
