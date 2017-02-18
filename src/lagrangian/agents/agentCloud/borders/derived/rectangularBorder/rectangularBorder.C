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

#include "rectangularBorder.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(rectangularBorder, 0);

addToRunTimeSelectionTable(borderModel, rectangularBorder, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
rectangularBorder::rectangularBorder
(
    Time& time,
    agentCloud& cloud,
    const dictionary& dict
)
:
    borderModel(time, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
//     wallModel()
{
    borderList_ = List<vectorList>(propsDict_.lookup("bordersList"));
    
    Info << "borderList = " << borderList_ << endl;
    
    // temporary
    treshold_ = 0.001;
    
    checkClosedEndedBorders();
    
    scalingValue_ = 20;
    
    if (propsDict_.found("scalingValue"))
    {
        scalingValue_ = readLabel(propsDict_.lookup("scalingValue"));
    }    

    // add wall force
    addWallForce_ = false;
    
    if (propsDict_.found("addWallForce"))
    {
        addWallForce_ = Switch(propsDict_.lookup("addWallForce"));
    }

    if(addWallForce_)
    {
        wallModel_ = autoPtr<agentWallForce>
        (
            agentWallForce::New(time, propsDict_)
        );
    }
   
   
    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );
    
    agentIds_ = ids.agentIds();
    
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rectangularBorder::~rectangularBorder()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void rectangularBorder::initialConfiguration()
{
    // initialise borders
    initialiseBorders();
    
    // delete agents inside borders
    
    const List< DynamicList<agent*> >& cellOccupancy
        = cloud_.cellOccupancy();
    
    DynamicList<agent*> agentsToDel;
    
    forAll(cellOccupancy, c)
    {
        const List<agent*>& molsInCell = cellOccupancy[c];

        forAll(molsInCell, mIC)
        {
            agent* molI = molsInCell[mIC];
            
            if(!molI->frozen())
            {
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
    }    
    
    forAll(agentsToDel, i)
    {
        cloud_.removeMolFromCellOccupancy(agentsToDel[i]);
        cloud_.deleteParticle(*agentsToDel[i]);
    }
}

void rectangularBorder::afterMove()
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
                if(findIndex(agentIds_, molI->id()) != -1)
                {
                    label index = interactionList_[c][i];
                    
                    if(isPointWithinBorder(index, molI->position()))
                    {
                        reflect(index, molI);
                    }
                }
            }
        }
    }    
    
    // rebuild cell Occupancy 
    cloud_.buildCellOccupancy();
}

void rectangularBorder::afterForce()
{
    if(addWallForce_)
    {
        const List< DynamicList<agent*> >& cellOccupancy
            = cloud_.cellOccupancy();

        forAll(cellOccupancy, c)
        {
            const List<agent*>& agentsInCell = cellOccupancy[c];

            forAll(agentsInCell, a)
            {
                agent* agentI = agentsInCell[a];
                
                if(findIndex(agentIds_, agentI->id()) != -1)
                {   
                    vector rI = agentI->position();
                    
//                     label closestCorner1 = -1;
//                     label closestBorder = -1;
//                     scalar dIWMin = GREAT;
                    
                    forAll(interactionList_[c], i)
                    {
                        label index = interactionList_[c][i];
                        
                        // find the closest edge to interact with                
                        for (int corner = 0; corner < borderList_[index].size()-1; corner++)
                        {
                            const vector& v1 = borderList_[index][corner];             
                            const vector& v2 = borderList_[index][corner+1];
                            
                            scalar edgeLength = mag(v2-v1);
                            vector v21 = v2-v1;
                            
                            // what's t?
                            scalar t = max(0, min(1, ((rI-v1) & v21)/(edgeLength*edgeLength)));        
                            
                            vector closestPointOnBorder = v1 + t*v21;
                            
                            // apply force 
                            agentI->f() += wallModel_->force(agentI, closestPointOnBorder);
                        }
                    }
                }
            }
        }    
    }
}


void rectangularBorder::checkClosedEndedBorders()
{
    forAll(borderList_, i)
    {
        const vector& rI = borderList_[i][0];
        const vector& rJ = borderList_[i][borderList_[i].size()-1];
        scalar rD = mag(rI - rJ);
        
        if(rD > treshold_)
        {
            FatalErrorIn("rectangularBorder::checkClosedEndedBorders()") 
                << "    This border has its end points not connected: " << nl
                << borderList_[i] << nl
                << "    Make the start point and end points the same." 
                << nl << abort(FatalError);  
        }
    }
}

void rectangularBorder::initialiseBorders()
{
    interactionList_.setSize(mesh_.nCells());
    
    //changeZ
    boundedBox bbMesh( mesh_.bounds().min(), mesh_.bounds().max());
    scalar Zmin=bbMesh.min().z();
    scalar Zmax=bbMesh.max().z();
    
    scalar dZ=mag(Zmax-Zmin)*0.1;
    
    Zmin -= dZ;
    Zmax += dZ;
    
    label nBorders = borderList_.size();
    
    normalVectors_.setSize(nBorders);
    midPoints_.setSize(nBorders);
    
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

    }
        
    Info << "midpoints = "<<  midPoints_ << endl;
    Info << "normal = " << normalVectors_ << endl;
    
//     Info << "interactionList = " << interactionList_ << endl;
}

bool rectangularBorder::isPointWithinBorder(label index, const vector& r)
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

// METHOD 1 
void rectangularBorder::reflect(label index, agent* p)
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
            FatalErrorIn("rectangularBorder::reflect()") 
                << "    No edge found = " << p->position()
                << "; What happened here?"
                << nl << abort(FatalError);   
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
        
        vector oldPosition = p->position();
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
            FatalErrorIn("rectangularBorder::reflect()") 
                << "    molecule left mesh at position = " << p->position()
                << ", starting from position = " << oldPosition
                << "; Implement more robust algorithm!" 
                << nl << abort(FatalError);  
        }
    }
}

// void rectangularBorder::reflect(label index, agent* p)
// {
//     // back track molecule 
//     scalar deltaT = mesh_.time().deltaT().value();
//     
//     vector r1 = p->position() - p->v()*deltaT;
//     vector r2 = p->position();
//     
//     Info << "r1 = " << r1 << endl;
//     Info << "r2 = " << r2 << endl;
//     
//     scalar x1 = r1.x();
//     scalar x2 = r2.x();
//     scalar y1 = r1.y();
//     scalar y2 = r2.y();
//     
//     vector n12 = r2 - r1;
//     scalar r12Mag = mag(n12);
//     n12 /= r12Mag;
//     
//     DynamicList<vector> midPoints;
//     DynamicList<vector> normals;    
//     DynamicList<vector> unitVectors;    
//     // loop over all line segments of border
//     
//     for (int i = 0; i < borderList_[index].size()-1 ; i++)
//     {
//         vector r3 = borderList_[index][i]; // vector of point 1
//         vector r4 = borderList_[index][i+1]; // vector of point 2
// 
//         Info << "r3 = " << r3 << endl;
//         Info << "r4 = " << r4 << endl;
//         
//         scalar x3 = r3.x();
//         scalar x4 = r4.x();
//         scalar y3 = r3.y();
//         scalar y4 = r4.y();
//         
//         scalar Px = ( (x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4) ) / 
//                     ( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) );
//                     
//         scalar Py = ( (x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4) ) /
//                     ( (x1-x2)*(y3-y4) - (y1 - y2)*(x3-x4)  );
//         
//         // midpoint (from Wikipedia)
//         vector rM = vector(Px, Py, 0.0);
//                 
//         Info << "rM = " << rM << endl;
//         
//         // test if the midpoint lies on both lines
//         
//         scalar rD1 = n12 & (rM - r1);
//         
//         Info << "rD1 = " << rD1 << endl;
//         
//         if( (rD1 >= 0) && (rD1 <= r12Mag) )
//         {
//             vector n34 = r4 - r3;
//             scalar r34Mag = mag(n34);
//             n34 /= r34Mag;
//             
//             scalar rD3 = n34 & (rM - r3);
//             Info << "rD3 = " << rD3 << endl;    
//             
//             if( (rD3 >= 0) && (rD3 <= r34Mag) )
//             {
//                 // found intersecting point
//                 unitVectors.append(n34);
//                 midPoints.append(rM);
//                 normals.append(normalVectors_[index][i]);
//             }
//         }
//     }
//     
//     Info << "midPoints = " <<  midPoints << endl;
//     
//     vector midPoint = vector::zero;
//     vector normal = vector::zero;
//     vector unitVector = vector::zero;
//     
//     if(midPoints.size() > 0)
//     {
//         scalar rDMax = GREAT;
//         
//         forAll(midPoints, i)
//         {
//             vector rM = midPoints[i];
//             scalar rD1 = n12 & (rM - r1);
//             
//             if(rD1 < rDMax)
//             {
//                 midPoint = rM;
//                 normal = normals[i];
//                 unitVector = unitVectors[i];
//                 rDMax = rD1;
//             }
//         }        
//     }
//     else
//     {
//         FatalErrorIn("rectangularBorder::reflect()") 
//                 << "    Algorithm Failure for agent at position = " << p->position()
//                 << nl << abort(FatalError);  
//     }
//     
//  
//     // return position to midPoint - offset 
//     
//     vector offset = n12*r12Mag/20;
//     
//     // new position:
//     
//     Info << "old position = " << p->position() << endl;
//     
//     p->position() = midPoint - offset;
//     
//     Info << "new position = " << p->position() << endl;
//     
//     // new velocity:
//     Info << "old vel = " << p->v() << endl;
//     
//     p->v() = (p->v() & unitVector)*unitVector;
//     
//     Info << "new vel = " << p->v() << endl;
//     
//     // check that you moved cell 
//     label cell = -1;
//     label tetFace = -1;
//     label tetPt = -1;
// 
//     mesh_.findCellFacePt
//     (
//         p->position(),
//         cell,
//         tetFace,
//         tetPt
//     );
//     
//     if(cell != -1)
//     {
//         if(p->cell() != cell)
//         {
//             p->cell() = cell;
// //                 Info << "changed cell to = " << cell << endl;        
//         }
//     }
//     else
//     {
//         FatalErrorIn("rectangularBorder::reflect()") 
//             << "    molecule left mesh at position = " << p->position()
//             << ", starting from position = " << r1
//             << "; Implement more robust algorithm!" 
//             << nl << abort(FatalError);  
//     }
// 
// }

// METHOD 3
// void rectangularBorder::reflect(label index, agent* p)
// {
//     // back track molecule 
//     scalar deltaT = mesh_.time().deltaT().value();
//     
//     vector r1 = p->position() - p->v()*deltaT;
//     vector r2 = p->position();
//     vector n12 = r2 - r1;
//     scalar r12Mag = mag(n12);
//     n12 /= r12Mag;    
//     
// 
//     
//     scalar rDMin = GREAT;
// //     label edge = -1;
//     vector u = vector::zero;
//     
//     forAll(midPoints_[index], i)
//     {
//         const vector& m = midPoints_[index][i];
//         const vector& n = normalVectors_[index][i];
//         
//         vector unitVector = borderList_[index][i]-borderList_[index][i+1];
//         unitVector /= mag(unitVector);
//         
//         scalar rD = mag((p->position()-m) & n);
//         
//         if( (rD < rDMin) )
//         {
//             rDMin = rD;
// //             edge = i;
//             u = unitVector;
//         }
//     }
//     
//     label value = scalingValue_; 
//     
//     vector r = r2;    
//     scalar dR = r12Mag/scalar(value);
//     
//     for (int i = 0; i < value+1 ; i++)
//     {
//         if(isPointWithinBorder(index, r))
//         {     
//             r -= i*dR*n12;            
//         }
//         else
//         {
//             break;
//         }
//     }
//     
//     // new position:
//     
// //     Info << "old position = " << p->position() << endl;
//     
//     p->position() = r;
//     
// //     Info << "new position = " << p->position() << endl;
//     
//     // new velocity:
// //     Info << "old vel = " << p->v() << endl;
//     
//     p->v() = mag(p->v())*u;
//     
// //     Info << "new vel = " << p->v() << endl;
//     
//     // check that you moved cell 
//     label cell = -1;
//     label tetFace = -1;
//     label tetPt = -1;
// 
//     mesh_.findCellFacePt
//     (
//         p->position(),
//         cell,
//         tetFace,
//         tetPt
//     );
//     
//     if(cell != -1)
//     {
//         if(p->cell() != cell)
//         {
//             p->cell() = cell;
// //                 Info << "changed cell to = " << cell << endl;        
//         }
//     }
//     else
//     {
//         FatalErrorIn("rectangularBorder::reflect()") 
//             << "    molecule left mesh at position = " << p->position()
//             << ", starting from position = " << r1
//             << "; Implement more robust algorithm!" 
//             << nl << abort(FatalError);  
//     }
// 
// }

void rectangularBorder::write(const fileName& pathName)
{
 
}

} // End namespace Foam

// ************************************************************************* //
