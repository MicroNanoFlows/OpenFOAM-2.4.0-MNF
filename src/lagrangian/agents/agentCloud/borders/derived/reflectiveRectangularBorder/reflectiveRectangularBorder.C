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

#include "reflectiveRectangularBorder.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(reflectiveRectangularBorder, 0);

addToRunTimeSelectionTable(borderModel, reflectiveRectangularBorder, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
reflectiveRectangularBorder::reflectiveRectangularBorder
(
    Time& time,
    agentCloud& cloud,
    const dictionary& dict
)
:
    borderModel(time, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
{
    borderList_ = List<vectorList>(propsDict_.lookup("bordersList"));
    
    Info << "borderList = " << borderList_ << endl;
    
    // temporary
    treshold_ = 0.001;
    
    checkClosedEndedBorders();

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

reflectiveRectangularBorder::~reflectiveRectangularBorder()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void reflectiveRectangularBorder::initialConfiguration()
{
    // initialise borders
    initialiseBorders();
    
    // delete molecules inside borders
    
    const List< DynamicList<agent*> >& cellOccupancy
        = cloud_.cellOccupancy();
    
    DynamicList<agent*> agentsToDel;
    
    forAll(cellOccupancy, c)
    {
        const List<agent*>& molsInCell = cellOccupancy[c];

        forAll(molsInCell, mIC)
        {
            agent* molI = molsInCell[mIC];
            
            forAll(interactionList_[c], i)
            {
                label index = interactionList_[c][i];
                
                if(isPointWithinBorder(index, molI->position()))
                {
                   agentsToDel.append(molI);
                }
            }
        }
    }    
    
    forAll(agentsToDel, i)
    {
        cloud_.removeMolFromCellOccupancy(agentsToDel[i]);
        cloud_.deleteParticle(*agentsToDel[i]);
    }
}

void reflectiveRectangularBorder::afterMove()
{
   // e.g. check for particles within borders
    
    const List< DynamicList<agent*> >& cellOccupancy
        = cloud_.cellOccupancy();

    forAll(cellOccupancy, c)
    {
        const List<agent*>& molsInCell = cellOccupancy[c];

        forAll(molsInCell, mIC)
        {
            agent* molI = molsInCell[mIC];
            
            forAll(interactionList_[c], i)
            {
                
                label index = interactionList_[c][i];
                
                if(isPointWithinBorder(index, molI->position()))
                {
//                     Info << "Warning: inside = " << molI->position() << endl;
                    
                    reflect(index, molI);
                }
            }
        }
    }    
}

void reflectiveRectangularBorder::afterForce()
{
/*    const List< DynamicList<agent*> >& cellOccupancy
        = cloud_.cellOccupancy();

    forAll(cellOccupancy, c)
    {
        const List<agent*>& molsInCell = cellOccupancy[c];

        forAll(molsInCell, mIC)
        {
            agent* molI = molsInCell[mIC];
            
            forAll(interactionList_[c], i)
            {
                label index = interactionList_[c][i];
                
                force 
                
            }
        }
    }*/       
}

void reflectiveRectangularBorder::checkClosedEndedBorders()
{
    forAll(borderList_, i)
    {
        const vector& rI = borderList_[i][0];
        const vector& rJ = borderList_[i][borderList_[i].size()-1];
        scalar rD = mag(rI - rJ);
        
        if(rD > treshold_)
        {
            FatalError
                << "reflectiveRectangularBorder::checkClosedEndedBorders() " << endl
                << "    This border has its end points not connected: " << nl
                << borderList_[i] << nl
                << "    Make the start point and end points the same." << endl;
        }
    }
}

void reflectiveRectangularBorder::initialiseBorders()
{
    interactionList_.setSize(mesh_.nCells());
    
    //changeZ
    boundedBox bbMesh( mesh_.bounds().min(), mesh_.bounds().max());
    scalar Zmin=bbMesh.min().z();
    scalar Zmax=bbMesh.max().z();
    
    scalar dZ=mag(Zmax-Zmin)*0.1;
    
    Zmin -= dZ;
    Zmax += dZ;
    
    normalVectors_.setSize(borderList_.size());
    midPoints_.setSize(borderList_.size());
    
    forAll(borderList_, i)
    {
        const vectorList& borderI = borderList_[i];
        
        vector min=vector(GREAT,GREAT,Zmin);
        vector max=vector(0.0,0.0,Zmax);
        
        // produce bounded box to border
        forAll(borderI, j)
        {
            const vector& vI = borderI[j];
            
            if(vI.x() < min.x())
            {
                min.x() = vI.x();
            }
            if(vI.y() < min.y())
            {
                min.y() = vI.y();
            }        
       
            if(vI.x() > max.x())
            {
                max.x() = vI.x();
            }        
            if(vI.y() > max.y())
            {
                max.y() = vI.y();
            }        
        }

        boundedBox bb (min,max);
        
        // scale box
        scalar rCut = cloud_.f().rCut();
        
        bb.expandII(rCut);
        
        // select cells 
        for (int c = 0; c < mesh_.nCells(); c++)
        {   
            const labelList& points = mesh_.cellPoints()[c];
            const labelList& faces = mesh_.cells()[c];

            label sizeOfList = points.size() + faces.size();

            pointField vectorPoints(sizeOfList, vector::zero);

            label counter = 0;

            forAll(points, p)
            {
                vectorPoints[counter] = mesh_.points()[points[p]];
                counter++;
            }
            
            forAll(faces, f)
            {
                vectorPoints[counter] = mesh_.faceCentres()[faces[f]];
                counter++;
            }
            
            bool inside = false;
            
            forAll(vectorPoints, v)
            {
                if(bb.contains(vectorPoints[v]))
                {
                    inside = true;
                }
            }
            
            if(inside)
            {
                interactionList_[c].append(i);
            }
        }
        
        label sz=borderList_[i].size()-1;
        normalVectors_[i].setSize(sz);
        midPoints_[i].setSize(sz);
        vector centrePoint = vector::zero;
        
        // set midpoints
        for (int p = 0; p < sz; p++)
        {
            const vector& vI=borderList_[i][p];
            centrePoint += vI;
        }
        
        centrePoint /= scalar(sz);
        
        Info << "centrePoint = " << centrePoint << endl;
        
        // set normals
        for (int p = 0; p < sz; p++)
        {
            const vector& vI=borderList_[i][p];
            const vector& vJ=borderList_[i][p+1];
            vector vIJ=vJ-vI;
            vector midPoint = vI + vIJ*0.5;
            midPoints_[i][p]=midPoint;
            
            vector n = midPoint-centrePoint; // pointing outwards
            n /= mag(n);
            
            // deal with skewness
            scalar theta = acos(vIJ & n / mag(vIJ));

            n *= sin(theta);
            n /= mag(n);
            
            normalVectors_[i][p] = n;
        }
        
        Info << "midpoints = "<<  midPoints_ << endl;
        Info << "normal = " << normalVectors_ << endl;
    }
    
//     Info << "interactionList = " << interactionList_ << endl;
}

bool reflectiveRectangularBorder::isPointWithinBorder(label index, const vector& r)
{
    bool inside = true;

    // a simple algorithm is implemented here
    // applicable only to rectangles and shapes like spheres or triangles
    
    forAll(midPoints_[index], i)
    {
        scalar rD = (r - midPoints_[index][i]) & normalVectors_[index][i];
        
        if(rD <= 0)
        {
           inside *= true;
        }
        else
        {
            inside *= false;
        }
    }
    
    return inside;
}

void reflectiveRectangularBorder::reflect(label index, agent* p)
{
    // find closest edge with normal in opposite direction
    
    scalar rDMin = GREAT;
    label edge = -1;
    
    forAll(midPoints_[index], i)
    {
        const vector& m = midPoints_[index][i];
        const vector& n = normalVectors_[index][i];
        
        scalar rD = mag((p->position()-m) & n);
        
        if( (rD < rDMin) && ( (p->v() & n) <= 0) )
        {
            rDMin = rD;
            edge = i;
        }
    }
    
    
    if (edge == -1)
    {
            FatalError
                << "reflectiveRectangularBorder::reflect() " << endl
                << "    No edge found = " << p->position()
                << "; What happened here?" << endl;  
    }
    else
    {
//         Info << "before velocity = " << p->v() << endl;
        
        const vector& n = normalVectors_[index][edge];
//         Info << "midpoint = " << midPoints_[index][edge] << endl;
//         Info << "normal = " << n << endl;
        
        scalar vN = ( p->v() & -n );
        
//         Info << "vN = " << vN << endl;
        
        p->v() += 2.0*vN*n;
//         Info << "after velocity = " << p->v() << endl;

          
        // move 
        
    
//         Info << "before position = " << p->position() << endl;        
        p->position() += 2.0*rDMin*n;
//         Info << "after position = " << p->position() << endl;        
        
        // check that you moved cell 
        label cell = -1;
        label tetFace = -1;
        label tetPt = -1;

        mesh_.findCellFacePt
        (
            p->position(),
            cell,
            tetFace,
            tetPt
        );
        
        if(cell != -1)
        {
            if(p->cell() != cell)
            {
                p->cell() = cell;
//                 Info << "changed cell to = " << cell << endl;        
            }
        }
        else
        {
            FatalError
                << "reflectiveRectangularBorder::reflect() " << endl
                << "    molecule left mesh at position = " << p->position()
                << "; Implement more robust algorithm!" << endl;            
        }
    }
}

void reflectiveRectangularBorder::write(const fileName& pathName)
{
 
}

} // End namespace Foam

// ************************************************************************* //
