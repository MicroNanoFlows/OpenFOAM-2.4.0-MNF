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

#include "hogmanayControl_Sol1.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(hogmanayControl_Sol1, 0);

addToRunTimeSelectionTable(agentController, hogmanayControl_Sol1, dictionary);


void hogmanayControl_Sol1::setBoundBox
(
    const dictionary& propsDict,
    boundedBox& bb,
    const word& name 
)
{
    const dictionary& dict(propsDict.subDict(name));
    
    vector startPoint = dict.lookup("startPoint");
    vector endPoint = dict.lookup("endPoint");

    bb.resetBoundedBox(startPoint, endPoint);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
hogmanayControl_Sol1::hogmanayControl_Sol1
(
    Time& t,
    agentCloud& cloud,
    const dictionary& dict
)
:
    agentController(t,  cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
//     agentId_()
{
    writeInTimeDir_ = false;
    writeInCase_ = false;

    
    setBoundBox(propsDict_, boxBottom_, "boxBottom");
    setBoundBox(propsDict_, boxTop_, "boxTop");
    
    {
        topAgentIds_.clear();

        selectAgentIds ids
        (
            cloud_.cP(),
            propsDict_,
            "topAgentIds"
        );

        topAgentIds_ = ids.agentIds();
    }    
    
    {
        bottomAgentIds_.clear();

        selectAgentIds ids
        (
            cloud_.cP(),
            propsDict_,
            "bottomAgentIds"
        );

        bottomAgentIds_ = ids.agentIds();
    }   
    
    {
        rightAgentIds_.clear();

        selectAgentIds ids
        (
            cloud_.cP(),
            propsDict_,
            "rightAgentIds"
        );

        rightAgentIds_ = ids.agentIds();
    }    
    
    {
        leftAgentIds_.clear();

        selectAgentIds ids
        (
            cloud_.cP(),
            propsDict_,
            "leftAgentIds"
        );

        leftAgentIds_ = ids.agentIds();
    }    
    
    tau_ = readScalar(propsDict_.lookup("tau"));
    
    frac1_ = 0.1;
    frac1_ = readScalar(propsDict_.lookup("frac1"));
    frac2_ = readScalar(propsDict_.lookup("frac2"));
    frac3_ = readScalar(propsDict_.lookup("frac3"));
        
        // borders 
    
    borderList_ = List<vectorList>(propsDict_.lookup("bordersList"));
    
    dictionary dict1 = propsDict_.subDict("attractiveForceModel");
    
    attractiveForceModel_ = autoPtr<agentWallForce>
    (
        agentWallForce::New(t, dict1)
    ); 

    
    initialiseBorders();
    
    
//     X1_ = readScalar(propsDict_.lookup("X1"));    
//     X2_ = readScalar(propsDict_.lookup("X2")); 
//     X3_ = readScalar(propsDict_.lookup("X3")); 
//     
    
    Y_ = readScalar(propsDict_.lookup("Y"));
    YMax_ = readScalar(propsDict_.lookup("YMax"));
    XRight_ = readScalar(propsDict_.lookup("XRight"));
    XLeft_ = readScalar(propsDict_.lookup("XLeft"));
    
    fireworkTime_ = readScalar(propsDict_.lookup("fireworkTime"));    
    evacuateTime_ = readScalar(propsDict_.lookup("evacuateTime"));
    
    meanDesSpeedE_ = readScalar(propsDict_.lookup("meanDesiredSpeedEvac"));
    desSpeedRangeE_ = readScalar(propsDict_.lookup("desiredSpeedRangeEvac"));


    meanDesSpeedF_ = readScalar(propsDict_.lookup("meanDesiredSpeedFire"));
    desSpeedRangeF_ = readScalar(propsDict_.lookup("desiredSpeedRangeFire"));    
    
//     deltaT_ = time_.deltaT().value();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

hogmanayControl_Sol1::~hogmanayControl_Sol1()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void hogmanayControl_Sol1::initialConfiguration()
{}

void hogmanayControl_Sol1::controlBeforeVelocityI()
{}

void hogmanayControl_Sol1::controlBeforeMove()
{}

void hogmanayControl_Sol1::controlBeforeForces()
{

}

void hogmanayControl_Sol1::controlDuringForces
(
    agent* molI,
    agent* molJ
)
{}

void hogmanayControl_Sol1::controlAfterForces()
{


    {
        IDLList<agent>::iterator p(cloud_.begin());
        
        for (p = cloud_.begin(); p != cloud_.end(); ++p)
        {
            if(p().eventTracker() == 0)
            {
                if(findIndex(topAgentIds_, p().id()) != -1)
                {
                    if(!boxTop_.contains(p().position()))
                    {
                        vector desDir = vector(0.0, (boxTop_.midpoint() - p().position()).y(), 0.0);
                        vector force = (p().desiredSpeed()*(desDir/mag(desDir)) - p().v())*p().mass() / tau_; 
                        p().f() += force;
/*                        vector R = vector(cloud_.rndGen().scalar01(), cloud_.rndGen().scalar01(), 0.0);
                        R /= mag(R);
                        p().f() += R*frac1_*mag(force); */        
                    }
                    else
                    {
//                         vector force = (p().desiredSpeed()*desDir - p().v())*p().mass() / tau_; 
//                         p().f() += force;
/*                        vector R = vector(cloud_.rndGen().scalar01(), cloud_.rndGen().scalar01(), 0.0);
                        R /= mag(R);
                        p().f() += R*frac2_; */                         
                    }
                }
                
                if(findIndex(bottomAgentIds_, p().id()) != -1)
                {
                    if(!boxBottom_.contains(p().position()))
                    {
                        vector desDir = vector::zero;
                        
                        if(p().position().y() > Y_)
                        {
                            desDir = vector(0.0, (boxBottom_.midpoint() - p().position()).y(), 0.0);
                        }
                        else 
                        {
                            desDir = vector((boxBottom_.midpoint() - p().position()).x(), 0.0, 0.0);
                        }
                        
                        vector force = (p().desiredSpeed()*(desDir/mag(desDir)) - p().v())*p().mass() / tau_; 
                        p().f() += force;
/*                        vector R = vector(cloud_.rndGen().scalar01(), cloud_.rndGen().scalar01(), 0.0);
                        R /= mag(R);
                        p().f() += R*frac1_*mag(force);  */       
                    }
                    else
                    {
//                         vector force = (p().desiredSpeed()*desDir - p().v())*p().mass() / tau_; 
//                         p().f() += force;
/*                        vector R = vector(cloud_.rndGen().scalar01(), cloud_.rndGen().scalar01(), 0.0);
                        R /= mag(R);
                        p().f() += R*frac2_; */                         
                    }                
                }
            }
            
            
            if( (time_.timeOutputValue() > fireworkTime_) && (p().eventTracker() == 0) )
            {
                p().eventTracker() = 1;
                p().desiredSpeed() = gaussianDistribution(meanDesSpeedF_, desSpeedRangeF_);
            }
            
            if(p().eventTracker() == 1)
            {
                if(!boxBottom_.contains(p().position()))
                {
                    vector desDir = vector::zero;
                    
                    if(p().position().y() > Y_ && p().position().y() < YMax_)
                    {
                        desDir = vector(0.0, (boxBottom_.midpoint() - p().position()).y(), 0.0);
                    }
                    else if(p().position().y() > YMax_ && p().position().x() < XRight_ && p().position().x() > XLeft_)
                    {
                        desDir = vector(0.0, (boxBottom_.midpoint() - p().position()).y(), 0.0);
                    }
                    else if(p().position().y() > YMax_ && p().position().x() > XRight_)
                    {
                        desDir = vector((boxBottom_.midpoint() - p().position()).x(), 0.0, 0.0);
                    }
                    else if(p().position().y() > YMax_ && p().position().x() < XLeft_)
                    {
                        desDir = vector((boxBottom_.midpoint() - p().position()).x(), 0.0, 0.0);
                    }
                    else 
                    {
                        desDir = vector((boxBottom_.midpoint() - p().position()).x(), 0.0, 0.0);
                    }
                        
                    vector force = (p().desiredSpeed()*(desDir/mag(desDir)) - p().v())*p().mass() / tau_; 
                    p().f() += force;
/*                    vector R = vector(cloud_.rndGen().scalar01(), cloud_.rndGen().scalar01(), 0.0);
                    R /= mag(R);
                    p().f() += R*frac1_*mag(force)*/;         
                }
                else
                {
/*                    vector R = vector(cloud_.rndGen().scalar01(), cloud_.rndGen().scalar01(), 0.0);
                    R /= mag(R);
                    p().f() += R*frac2_;   */                       
                }
            }
            
            if((time_.timeOutputValue() > evacuateTime_) && (p().eventTracker() == 1))
            {
                p().eventTracker() = 2;
                p().desiredSpeed() = gaussianDistribution(meanDesSpeedE_, desSpeedRangeE_);
            }    
            
            
            // evacuation
            if(p().eventTracker() == 2)            
            {             
                
                vector desDir = vector::zero;                    
                
                if(p().position().y() > Y_ && p().position().y() < YMax_)
                {
                    desDir = vector(0.0, 1.0, 0.0);
                }
                else if(p().position().y() < Y_ && p().position().y() > YMin2_) 
                {
                    desDir = vector(1.0, 0.0, 0.0);
                }
                else if(p().position().y() < YMin1_) 
                {
                    desDir = vector(-1.0, 0.0, 0.0);
                }
                else if(p().position().y() < YMin_) 
                {
                    desDir = vector(-1.0, 0.0, 0.0);
                }
                    
                vector force = (p().desiredSpeed()*(desDir/mag(desDir)) - p().v())*p().mass() / tau_; 
                p().f() += force;
                vector R = vector(cloud_.rndGen().scalar01(), cloud_.rndGen().scalar01(), 0.0);
                R /= mag(R);
                p().f() += R*frac3_*mag(force);         
            }
        }
    }
    
    
    // forces from borders
    {
        const List< DynamicList<agent*> >& cellOccupancy
            = cloud_.cellOccupancy();

        forAll(cellOccupancy, c)
        {
            const List<agent*>& agentsInCell = cellOccupancy[c];

            forAll(agentsInCell, a)
            {
                agent* agentI = agentsInCell[a];
                
                if(agentI->eventTracker() == 1)
                {                
                    vector rI = agentI->position();
                
                    scalar minD = GREAT;
                    vector pointOnBorder = vector::zero;
                    
                    label borderNumber = 0;
                    
                    // find the closest edge to interact with                
                    for (int corner = 0; corner < borderList_[borderNumber].size()-1; corner++)
                    {
                        const vector& v1 = borderList_[borderNumber][corner];             
                        const vector& v2 = borderList_[borderNumber][corner+1];
                        
                        scalar edgeLength = mag(v2-v1);
                        vector v21 = v2-v1;
                        
                        scalar t = max(0, min(1, ((rI-v1) & v21)/(edgeLength*edgeLength)));        
                        
                        vector closestPointOnBorder = v1 + t*v21;
                        vector rWI = closestPointOnBorder - rI;
                        scalar dWI = mag(rWI);
                    
                        if(dWI < minD)
                        {
                            minD = dWI;
                            pointOnBorder = closestPointOnBorder;
                        }                            
                    }
                    
                    // apply attractive force 
                    agentI->f() += attractiveForceModel_->force(agentI, pointOnBorder);                        

                }
            }
        }   
    }
}


void hogmanayControl_Sol1::initialiseBorders()
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
    
//     sideTop_.setSize(3, vector::zero);
//     sideBottom_.setSize(3, vector::zero);    
    
    centrePoints_.setSize(3, vector::zero);    
    
    forAll(borderList_, i)
    {
        const vectorList& borderI = borderList_[i];
        
        vector min=vector(GREAT,GREAT,Zmin);
        vector max=vector(0.0,0.0,Zmax);
        vector centrePoint = vector::zero;
        
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
            
            centrePoint += vI;
        }

        boundedBox bb (min,max);
        
            
        centrePoints_[i] = centrePoint/scalar(borderI.size());
        
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


//         label sz=borderList_[i].size()-1;
//         normalVectors_[i].setSize(sz);
//         midPoints_[i].setSize(sz);
//         vector centrePoint = vector::zero;
//         
//         // set midpoints
//         for (int p = 0; p < sz; p++)
//         {
//             const vector& vI=borderList_[i][p];
//             centrePoint += vI;
//         }
//         
//         centrePoint /= scalar(sz);
//         
//         Info << "centrePoint = " << centrePoint << endl;
//         
//         // set normals
//         for (int p = 0; p < sz; p++)
//         {
//             const vector& vI=borderList_[i][p];
//             const vector& vJ=borderList_[i][p+1];
//             vector vIJ=vJ-vI;
//             vector midPoint = vI + vIJ*0.5;
//             midPoints_[i][p]=midPoint;
//             
//             vector n = midPoint-centrePoint; // pointing outwards
//             n /= mag(n);
//             
//             // deal with skewness
//             scalar theta = acos(vIJ & n / mag(vIJ));
// 
//             n *= sin(theta);
//             n /= mag(n);
//             
//             normalVectors_[i][p] = n;
//         }

    }
        
//     Info << "midpoints = "<<  midPoints_ << endl;
//     Info << "normal = " << normalVectors_ << endl;    


}




void hogmanayControl_Sol1::controlAfterVelocityII()
{}

void hogmanayControl_Sol1::calculateProperties()
{}

void hogmanayControl_Sol1::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
    }
}

scalar hogmanayControl_Sol1::gaussianDistribution
(
    const scalar& mean,
    const scalar& range
 
)
{
    return mean + 0.25*range*cloud_.rndGen().GaussNormal();
}


} // End namespace Foam

// ************************************************************************* //
