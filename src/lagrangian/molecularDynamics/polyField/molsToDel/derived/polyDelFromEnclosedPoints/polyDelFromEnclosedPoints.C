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

#include "polyDelFromEnclosedPoints.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyDelFromEnclosedPoints, 0);

addToRunTimeSelectionTable(polyMolsToDeleteModel, polyDelFromEnclosedPoints, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyDelFromEnclosedPoints::polyDelFromEnclosedPoints
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyMolsToDeleteModel(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    molIds_()
{
    
    molIds_.clear();

    selectIds ids
    (
        molCloud_.pot(),
        propsDict_
    );

    molIds_ = ids.molIds();

    
    // delete outside list of points 
    
    deleteExternalMols_ = true;
 
    if (propsDict_.found("deleteOutside"))
    {
        deleteExternalMols_ = Switch(propsDict_.lookup("deleteOutside"));  
    }
    
    // read in list of points 
    // this has to be an ordered list - points need to be after each other
    List<vector> points = List<vector>(propsDict_.lookup("points"));
    
    vector translate = vector::zero;
    
    if (propsDict_.found("translate"))
    {    
        translate = propsDict_.lookup("translate");
    }
    
    normal_ = propsDict_.lookup("normal");
    
    points_.setSize(points.size(), vector::zero);
    
    forAll(points, i)
    {
        points_[i] = points[i]+translate;
    }
    Info << "points " << points_ << endl;
    
    if(!(points_.size() > 0))
    {
        FatalErrorIn("polyDelFromEnclosedPoints::polyDelFromEnclosedPoints()")
            << "Did not specify a list of points in dictionary!"
            << exit(FatalError);        
    }
    
    // Assume 2D for now and list is ordered
    
    // find centre point 
    
    centre_ = vector::zero;
    
    forAll(points_, i)
    {
        centre_ += points_[i];
    }
    
    centre_ /= scalar(points_.size());
    
    Info << "Centre = " << centre_ << endl;
    
    findMolsToDel();
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyDelFromEnclosedPoints::~polyDelFromEnclosedPoints()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool polyDelFromEnclosedPoints::isPointWithinTriangle
(
    const vector& r1,
    const vector& r2,
    const vector& rI
)
{
    bool inside = false;
    
    // face centres
    vector fC12=(r1+r2)/2;
    vector fC1C=(r1+centre_)/2;
    vector fC2C=(r2+centre_)/2;
    
    vector r1C = centre_ - r1;
    vector r12 = r2 - r1;
    vector r2C = centre_ - r2;
    
    vector fN12 = r12 ^ normal_;
    fN12 /= mag(fN12);
    
    vector fN1C = normal_ ^ r1C;
    fN1C /= mag(fN1C);
    
    vector fN2C = r2C ^ normal_;
    fN2C /= mag(fN2C);
    
    // solve for normals being inwards
    
    vector localC = (r1 + r2 + centre_)/3.0;
    
    if( ((localC-fC12) & fN12) < 0)
    {
        fN12 *= -1;
    }

    if( ((localC-fC1C) & fN1C) < 0)
    {
        fN1C *= -1;
    }

    if( ((localC-fC2C) & fN2C) < 0)
    {
        fN2C *= -1;
    }
    
//     Info << "triangle" << endl;
//     Info << "r1 = " << r1 << endl;
//     Info << "r2 = " << r2 << endl;
//     Info << "rI = " << rI << endl;
//     Info << "centre = " << centre_ << endl;
//     Info << "fC12 = " << fC12 << endl;
//     Info << "fC1C = " << fC1C << endl;
//     Info << "fC2C = " << fC2C << endl;
// 
//     Info << "fN12 = " << fN12 << endl;
//     Info << "fN1C = " << fN1C << endl;
//     Info << "fN2C = " << fN2C << endl;
    
    vector t1 = rI - fC12;
    vector t2 = rI - fC1C;
    vector t3 = rI - fC2C;
    
    if
    (
        ((t1 & fN12) >= 0) && 
        ((t2 & fN1C) >= 0) &&
        ((t3 & fN2C) >= 0)
    )
    {
        inside = true;
//         Info << "in" << endl;
    }
//     else
//     {
//         Info << "out" << endl;
//     }
//     vector v0 = x1;
//     vector v1 = centre_-v0;
//     vector v2 = x2-v0;
//     vector v = rI;
//     vector newV = v ^ v2;
//     Info << "cross = " << newV << endl;
//     scalar a = ( (v ^ v2) - (v0 ^ v2) )/ (v1 ^ v2);
//     scalar b = -( (v ^ v1) - (v0 ^ v1) ) /(v1 ^ v2);
    
//     if(( a > 0 ) && (b > 0) && ((a+b) < 1) )
//     scalar a = mag(newV);
//     if( a > 0)
//     {
//         inside = true;
//     }
    
    return inside;
}

bool polyDelFromEnclosedPoints::isPointWithinRegion(const vector& rI)
{
    
//     Info << nl << "Is point in region = " << rI << endl;
    
    bool inside = false;
    
    forAll(points_, i)
    {
        label I = i;
        label J = i+1;
        
        if(J > points_.size()-1)
        {
            J = 0;
        }
        
        if
        (
            isPointWithinTriangle
            (
                points_[I],
                points_[J],
                rI
            )
        )
        {
            inside = true;
        }
    }
    
//     if (inside)
//     {
//         Info << "in" << endl;
//     }
//     else
//     {
//         Info << "out" <<  endl;
//     }
    
    return inside;
}

void polyDelFromEnclosedPoints::findMolsToDel()
{
    DynamicList<polyMolecule*> molsToDel;

    label initialSize = molCloud_.size();
    
    {    
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        
        for
        (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
        )
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                polyMolecule* molI = &mol();
                
                if(deleteExternalMols_)
                {
                    if(!isPointWithinRegion(molI->position()))
                    {
                        molsToDel.append(molI);
                    }
                }
                else
                {
                    if(isPointWithinRegion(molI->position()))
                    {
                        molsToDel.append(molI);
                    }
                }
            }
        }
    }
    
    //molsToDel.shrink();

    forAll(molsToDel, m)
    {
        deleteMolFromMoleculeCloud(*molsToDel[m]);
    }

    label molsKept = initialSize - molsToDel.size();

    Info<< tab << " initial polyMolecules: " <<  initialSize 
        << ", polyMolecules kept: " <<  molsKept
        << ", polyMolecules removed: " << molsToDel.size() 
        << endl;


    // as a precaution: rebuild cell occupancy
    molCloud_.rebuildCellOccupancy();
    molCloud_.prepareInteractions();
}


} // End namespace Foam

// ************************************************************************* //
