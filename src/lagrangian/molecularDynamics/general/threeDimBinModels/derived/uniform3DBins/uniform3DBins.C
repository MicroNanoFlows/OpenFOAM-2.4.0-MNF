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
    uniform3DBins

Description

\*----------------------------------------------------------------------------*/

#include "uniform3DBins.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(uniform3DBins, 0);

addToRunTimeSelectionTable(threeDimBinModel, uniform3DBins, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
uniform3DBins::uniform3DBins
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    threeDimBinModel(mesh, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    startPoint_(propsDict_.lookup("startPoint")),
    endPoint_(propsDict_.lookup("endPoint")),
    unitVector_((endPoint_ - startPoint_)/mag(endPoint_ - startPoint_)),
    rSEMag_(mag(endPoint_ - startPoint_)),
    nBins_(propsDict_.lookup("nBins")),
    binWidth_(mag(endPoint_ - startPoint_)/(nBins_))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

uniform3DBins::~uniform3DBins()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// cellI is a dummy variable
List<label> uniform3DBins::isPointWithinBin
(
    const vector& rI,
    const label& cellI
)
{
    List<label> binNumbers(3, -1);
/*
    vector rSI = rI - startPoint_;
    scalar rD = rSI & unitVector_;
    vector n(label(rD/binWidth_[0]), label(rD/binWidth_[1]), label(rD/binWidth_[2]));
    
    if
    (
        (rD <= rSEMag_) && (rD >= 0.0)
    )
    {
      forAll(n, i)
      {
        if(n[i] == nBins_[i])
        {
            n[i]--;
        }

        if( (n[i] >= 0) && (n[i] < nBins_[i]) )
        {
            binNumber[i] = n[i];
        }
      }
    }
*/
    return binNumbers;
}

vectorField uniform3DBins::binPositionsX()
{}

scalarField uniform3DBins::binPositionsY()
{}

scalarField uniform3DBins::binPositionsZ()
{}

vectorField uniform3DBins::binPositionsXYZ()
{
    scalar totalBins = nBins_[0]*nBins_[1]*nBins_[2];
    vectorField positions(totalBins, vector::zero);
    label i, j, k;

    for (i = 0; i < nBins_[0]; i++)
    {
        for (j = 0; j < nBins_[1]; j++)
        {
            for (k = 0; k < nBins_[2]; k++)
            {
                scalar pos_x = 0.5*binWidth_[i] + scalar(i)*binWidth_[i];
                scalar pos_y = 0.5*binWidth_[j] + scalar(j)*binWidth_[j];
                scalar pos_z = 0.5*binWidth_[k] + scalar(k)*binWidth_[k];

                vector pos(pos_x, pos_y, pos_z);

                positions[i*j*k] = pos;
            }
        }
    }

    return positions;
}

void uniform3DBins::write
(
    const fileName& path,
    const word& name
)
{}

vectorField uniform3DBins::position()
{
    vectorField positions(nBins_[0], vector::zero);

    /*
    forAll(positions, i)
    {
        positions[i] = startPoint_ + (0.5 + scalar(i))*binWidth_*unitVector_;
    }
    */

    return positions;
}

List<label> uniform3DBins::nBins()
{
    List<label> nBins(3, -1);

    nBins[0] = nBins_[0];
    nBins[1] = nBins_[1];
    nBins[2] = nBins_[2];

    return nBins;
}

scalar uniform3DBins::binVolume(const label& n)
{
  /*
  scalar volume = area_*binWidth_;

  return volume;
  */

  return 0;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
