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

#include "jamaratControl.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(jamaratControl, 0);

addToRunTimeSelectionTable(agentController, jamaratControl, dictionary);


void jamaratControl::setBoundBox
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
jamaratControl::jamaratControl
(
    Time& t,
    agentCloud& cloud,
    const dictionary& dict
)
:
    agentController(t,  cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    agentId_()
{
    writeInTimeDir_ = false;
    writeInCase_ = false;

    
    // insertion 
    
    {
        const List<word>& idList(cloud_.cP().agentIds());
        const word agentId = propsDict_.lookup("agentId");
        agentId_ = findIndex(idList, agentId);
    
        if(agentId_ == -1)
        {
            FatalErrorIn("jamaratControl::jamaratControl()")
                << "Cannot find agentId: " << agentId << nl << "in: "
                << time_.time().system()/"controllersDict"
                << exit(FatalError);
        }
    }
    
    N_ = readLabel(propsDict_.lookup("noOfAgents_everyWriteInterval"));

    meanMass_ = readScalar(propsDict_.lookup("meanMass"));
    massRange_ = readScalar(propsDict_.lookup("massRange"));

    meanRadius_ = readScalar(propsDict_.lookup("meanRadius"));
    radiusRange_ = readScalar(propsDict_.lookup("radiusRange"));

    meanDesSpeed_ = readScalar(propsDict_.lookup("meanDesiredSpeed"));
    desSpeedRange_ = readScalar(propsDict_.lookup("desiredSpeedRange"));
    
//     setBoundBox(propsDict_, samplingBox_, "samplingBox");
    setBoundBox(propsDict_, controlBox_, "controlBox");
    
    
    
    // will force information
    
    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();    
    
    desiredDirection_ = propsDict_.lookup("desiredDirection");
    tau_ = readScalar(propsDict_.lookup("tau"));
    
//     frac_ = 0.1;
//     frac2_ = 0.05;
    fracP_ = readScalar(propsDict_.lookup("fracPillar"));
    fracR_ = readScalar(propsDict_.lookup("fracRoad"));
    fracV_ = readScalar(propsDict_.lookup("fracVel"));        
        // borders 
    
    borderList_ = List<vectorList>(propsDict_.lookup("bordersList"));
    
    dictionary dict1 = propsDict_.subDict("attractiveForceModel");
    
    attractiveForceModel_ = autoPtr<agentWallForce>
    (
        agentWallForce::New(t, dict1)
    ); 
    
    dictionary dict2 = propsDict_.subDict("repulsiveForceModel");
    
    repulsiveForceModel_ = autoPtr<agentWallForce>
    (
        agentWallForce::New(t, dict2)
    );    
    
//     setBoundBoxes();
    
    initialiseBorders();
    
    
    X1_ = readScalar(propsDict_.lookup("X1"));    
    X2_ = readScalar(propsDict_.lookup("X2")); 
    X3_ = readScalar(propsDict_.lookup("X3")); 
    
    searchTime_ = readScalar(propsDict_.lookup("searchTime"));    
    throwingTime_ = readScalar(propsDict_.lookup("throwingTime"));
    patienceTime_ = readScalar(propsDict_.lookup("patienceTime"));    
    maxThrowingDistance_ = readScalar(propsDict_.lookup("maxThrowingDistance"));
    
    
    deltaT_ = time_.deltaT().value();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

jamaratControl::~jamaratControl()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void jamaratControl::initialConfiguration()
{}

void jamaratControl::controlBeforeVelocityI()
{}

void jamaratControl::controlBeforeMove()
{}

void jamaratControl::controlBeforeForces()
{
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
                    
                    reflect(index, molI);
                }
            }
        }
    }    
    
    // rebuild cell Occupancy 
    cloud_.buildCellOccupancy();    
    
}


// METHOD 1 
void jamaratControl::reflect(label index, agent* p)
{
    scalar minD = GREAT;
    vector pointOnBorder = vector::zero;
//     label vertex = 0;
    
    // find the closest edge to interact with                
    for (int corner = 0; corner < borderList_[index].size()-1; corner++)
    {
        const vector& v1 = borderList_[index][corner];             
        const vector& v2 = borderList_[index][corner+1];
        
        scalar edgeLength = mag(v2-v1);
        vector v21 = v2-v1;
        
        scalar t = max(0, min(1, ((p->position()-v1) & v21)/(edgeLength*edgeLength)));        
        
        vector closestPointOnBorder = v1 + t*v21;
        
        vector rWI = closestPointOnBorder - p->position();
        scalar dWI = mag(rWI);
       
        if(dWI < minD)
        {
            minD = dWI;
            pointOnBorder = closestPointOnBorder;
//             vertex = corner;
        }
    }
    
/*        Info << " agent position  = " << p->position() 
             << " closestPoint on border  = " << closestPointOnBorder
             << " border points = " << v1 << " & " << v2
             <<  " centre point = " << rC 
             << " distance from wall = " << dWI
             << endl;*/    

    vector rC = centrePoints_[index];
    vector rIC = rC - p->position();
    vector rWC = rC - pointOnBorder;
    vector rWI = pointOnBorder - p->position();
    scalar dWI = mag(rWI);
    
    
    // reflect 
    if((mag(rIC)-p->radius()) < mag(rWC))
    {
/*        Info << " agent position  = " << p->position() 
             << " closestPoint on border  = " << pointOnBorder
             << " border points = " << borderList_[index][vertex] << " & " << borderList_[index][vertex+1]
             <<  " centre point = " << rC 
             << " distance from wall = " << dWI
             << endl; */ 
        vector nij = rWI/dWI;
        
        vector newPosition = p->position() + 2.0*rWI;
        

        vector tij = vector (-nij.y(), nij.x(), 0);
        vector vNew = mag(p->v())*tij;
        p->v() = vNew;
        
        // move 
        
        vector oldPosition = p->position();
//             Info << "before position = " << p->position() << endl;        
        p->position() = newPosition;
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
                p->tetFace() = tetFace;
                p->tetPt() = tetPt;
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

void jamaratControl::controlDuringForces
(
    agent* molI,
    agent* molJ
)
{}

void jamaratControl::panic
(
    agent* p,
    const vector& willForce
)
{
    scalar desSpeed = mag(p->desiredSpeed());
    
    if(mag(p->v()) < desSpeed*fracV_)
    {
        p->fraction() += deltaT_;
        vector Fpanic = (p->fraction()/patienceTime_)*willForce;
        p->f() += Fpanic;
//         Info << "panic force = " << Fpanic << " time = " << p->t() << endl;
    }
    else
    {
        p->fraction() = 0.0;
    }
}

void jamaratControl::controlAfterForces()
{
    if(time_.outputTime())
    {
//         Info << "jamaratControl: control" << endl;
        
        insertAgents(N_);
    }
    

    {
        IDLList<agent>::iterator p(cloud_.begin());
        
        for (p = cloud_.begin(); p != cloud_.end(); ++p)
        {
            // will force 
            
            if(p().eventTracker() == 0)
            {
                if(p().position().x() >= X1_)
                {
                    p().eventTracker() = 1;
                }
                
                if(findIndex(agentIds_, p().id()) != -1)
                {
                    vector willForce = (p().desiredSpeed()*desiredDirection_ - p().v())*p().mass() / tau_; 
                    p().f() += willForce;
                    vector R = vector(0.0, 2.0*cloud_.rndGen().scalar01()-1.0, 0.0);
                    R /= mag(R);
                    p().f() += R*fracR_*mag(willForce);
                }
            }
            
            // check in throwing box    
//             if(p().eventTracker() == 1)            
//             {
//                 if(boxes_[0].contains(p().position()))
//                 {
//                     p().eventTracker() = 2;
//                     p().t() = searchTime_;
//                 }
//             }
            
            // search for a nice position
            if(p().eventTracker() == 2)
            {
                if(p().t() > 0)
                {
                    p().t() -= deltaT_;
                    
                    if(p().t() < 0)
                    {
                        p().t() = 0.0;
                        p().eventTracker() = 3;
                        p().t() = throwingTime_;
                    }
                }
            }
            
            // throw stones
            if(p().eventTracker() == 3)
            {
                p().v() = vector::zero;
                
                if(p().t() > 0)
                {
                    p().t() -= deltaT_;
                    
                    if(p().t() < 0)
                    {
                        p().t() = 0.0;
                        p().eventTracker() = 4;
                    }
                }
            }
            
            
            //- iterate to 2nd pillar
            
            // switch on will force 
            if(p().eventTracker() == 4)
            {
                if(p().position().x() >= X2_)
                {
                    p().t() = 0.0;
                    p().eventTracker() = 5;
                }
                
                if(p().position().x() < P_[0])
                {

                    vector willForce = (p().desiredSpeed()*desiredDirection_ - p().v())*p().mass() / tau_; 
                    p().f() += willForce;
                    vector R = vector(0.0, 2.0*cloud_.rndGen().scalar01()-1.0, 0.0);
                    R /= mag(R);
                    p().f() += R*fracP_*mag(willForce);                  
                    panic(&p(), willForce);
                    
                }
                else
                {
                    vector willForce = (p().desiredSpeed()*desiredDirection_ - p().v())*p().mass() / tau_; 
                    p().f() += willForce;
                    vector R = vector(0.0, 2.0*cloud_.rndGen().scalar01()-1.0, 0.0);
                    R /= mag(R);
                    p().f() += R*fracR_*mag(willForce);                    
                }
            }
            
            // check in throwing box 
//             if(p().eventTracker() == 5)            
//             {
//                 if(boxes_[1].contains(p().position()))
//                 {
//                     p().eventTracker() = 6;
//                     p().t() = searchTime_;
//                 }
//             }
            
            // search for a nice position
            if(p().eventTracker() == 6)            
            {
                if(p().t() > 0)
                {
                    p().t() -= deltaT_;
                    
                    if(p().t() < 0)
                    {
                        p().t() = 0.0;
                        p().eventTracker() = 7;
                        p().t() = throwingTime_;
                    }
                }
            }
            
            // throw stones
            if(p().eventTracker() == 7)            
            {
                p().v() = vector::zero;
                
                if(p().t() > 0)
                {
                    p().t() -= deltaT_;
                    
                    if(p().t() < 0)
                    {
                        p().t() = 0.0;
                        p().eventTracker() = 8;
                    }
                }                 
            }            
            

            //- iterate to 3nd pillar
            
            // switch on will force 
            if(p().eventTracker() == 8)
            {
                if(p().position().x() >= X3_)
                {
                    p().eventTracker() = 9;
                }
                
                if(p().position().x() < P_[1])
                {
                    vector willForce = (p().desiredSpeed()*desiredDirection_ - p().v())*p().mass() / tau_; 
                    p().f() += willForce;
                    vector R = vector(0.0, 2.0*cloud_.rndGen().scalar01()-1.0, 0.0);
                    R /= mag(R);
                    p().f() += R*fracP_*mag(willForce);
                    panic(&p(), willForce);
                }
                else
                {
                    vector willForce = (p().desiredSpeed()*desiredDirection_ - p().v())*p().mass() / tau_; 
                    p().f() += willForce;
                    vector R = vector(0.0, 2.0*cloud_.rndGen().scalar01()-1.0, 0.0);
                    R /= mag(R);
                    p().f() += R*fracR_*mag(willForce);                    
                }
            }
            
            // check in throwing box 
//             if(p().eventTracker() == 9)            
//             {
//                 if(boxes_[2].contains(p().position()))
//                 {
//                     p().eventTracker() = 10;
//                     p().t() = searchTime_;
//                 }
//             }
            
            // search for a nice position
            if(p().eventTracker() == 10)            
            {
                if(p().t() > 0)
                {
                    p().t() -= deltaT_;
                    
                    if(p().t() < 0)
                    {
                        p().t() = 0.0;
                        p().eventTracker() = 11;
                        p().t() = throwingTime_;
                    }
                }
            }
            
            // throw stones
            if(p().eventTracker() == 11)            
            {
                p().v() = vector::zero;
                
                if(p().t() > 0)
                {
                    p().t() -= deltaT_;
                    
                    if(p().t() < 0)
                    {
                        p().t() = 0.0;
                        p().eventTracker() = 12;
                    }
                }                 
            }

            // switch on will force 
            if(p().eventTracker() == 12)
            {
                if(p().position().x() < P_[2])
                {
                    vector willForce = (p().desiredSpeed()*desiredDirection_ - p().v())*p().mass() / tau_; 
                    p().f() += willForce;
                    vector R = vector(0.0, 2.0*cloud_.rndGen().scalar01()-1.0, 0.0);
                    R /= mag(R);
                    p().f() += R*fracP_*mag(willForce);           
                    panic(&p(), willForce);
                }
                else
                {
                    vector willForce = (p().desiredSpeed()*desiredDirection_ - p().v())*p().mass() / tau_; 
                    p().f() += willForce;
                    vector R = vector(0.0, 2.0*cloud_.rndGen().scalar01()-1.0, 0.0);
                    R /= mag(R);
                    p().f() += R*fracR_*mag(willForce);                    
                }
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
                
                if(findIndex(agentIds_, agentI->id()) != -1)
                {   
                    vector rI = agentI->position();
                    
                    
                    // attractive force 
                    
                    label borderNumber = -1;
                    
                    if(agentI->eventTracker() >= 1)
                    {
                        borderNumber = 0;
                    }
                    
                    if(agentI->eventTracker() == 4)
                    {
                        borderNumber = -1;
                    }
                    
                    if(agentI->eventTracker() >= 5)
                    {
                        borderNumber = 1;
                    }

                    if(agentI->eventTracker() == 8)
                    {
                        borderNumber = -1;
                    }
                    
                    if(agentI->eventTracker() >= 9)
                    {
                        borderNumber = 2;
                    }
                    
                    if(agentI->eventTracker() >= 12)
                    {
                        borderNumber = -1;
                    }                    
                    
                    if(borderNumber != -1)
                    {
                        scalar minD = GREAT;
                        vector pointOnBorder = vector::zero;
                        
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

                        // check in throwing box
                        if(agentI->eventTracker() == 1)
                        {
                            if(minD <= maxThrowingDistance_)
                            {
                                agentI->eventTracker() = 2;
                                agentI->t() = searchTime_;
                            }
                        }
                        
                        if(agentI->eventTracker() == 5)
                        {
                            if(minD <= maxThrowingDistance_)
                            {
                                agentI->eventTracker() = 6;
                                agentI->t() = searchTime_;
                            }
                        }
                        
                        if(agentI->eventTracker() == 9)
                        {
                            if(minD <= maxThrowingDistance_)
                            {
                                agentI->eventTracker() = 10;
                                agentI->t() = searchTime_;
                            }
                        }
                    }
                    
                    // repulsive force 
                    
                    forAll(interactionList_[c], i)
                    {
                        label index = interactionList_[c][i];

                        scalar minD = GREAT;
                        vector pointOnBorder = vector::zero;
                        
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
                            vector rWI = closestPointOnBorder - rI;
                            scalar dWI = mag(rWI);
                        
                            if(dWI < minD)
                            {
                                minD = dWI;
                                pointOnBorder = closestPointOnBorder;
                            }
                        }
                        
                        // apply repulsive force
                        agentI->f() += repulsiveForceModel_->force(agentI, pointOnBorder);                        
                    }
                }
            }
        }   
    }
}

void jamaratControl::insertAgents(const label& N)
{
    boundedBox& bb = controlBox_;
    
    DynamicList<vector> positions;
    
    for (label i = 0; i < N; i++)
    {
        vector pos = vector
                     ( 
                        (bb.max().x() - bb.min().x())*cloud_.rndGen().scalar01() + bb.min().x(), 
                        (bb.max().y() - bb.min().y())*cloud_.rndGen().scalar01() + bb.min().y(),  
                        0.0
                     );
        
        positions.append(pos);
    }    

    //*** deal with parallel processing here ***//
    
    
    scalar massI = gaussianDistribution(meanMass_, massRange_);
    scalar radius = gaussianDistribution(meanRadius_, radiusRange_);
    scalar desiredSpeed = gaussianDistribution(meanDesSpeed_, desSpeedRange_);
    
    
    label nAgentsInserted = 0;
    
    // insert agents
    
    forAll(positions, i)
    {
        label cell = -1;
        label tetFace = -1;
        label tetPt = -1;

        mesh_.findCellFacePt
        (
            positions[i],
            cell,
            tetFace,
            tetPt
        );
        
        if(cell != -1)
        {
            cloud_.createAgent
            (
                positions[i],
                cell,
                tetFace,
                tetPt,     
                vector::zero, //v
                vector::zero,
                vector::zero,
                vector::zero,
                massI, // mass 
                radius,
                desiredSpeed, 
                0.0,
                GREAT,
                1.0,
                0.0,
                agentId_,
                0,
                cloud_.getTrackingNumber()
            );
            
            nAgentsInserted++;
        }
    }
    
    Info << "inserted = " << nAgentsInserted << endl;
}

void jamaratControl::initialiseBorders()
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
    P_.setSize(3, 0.0);    
    
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
//         sideTop_[i]=centrePoints_[i].midpoint()+vector(0,12,0);
//         sideBottom_[i]=centrePoints_[i].midpoint()-vector(0,12,0);
        
        P_[i] = max.x();
        
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

// void jamaratControl::setBoundBoxes()
// {
//     PtrList<entry> boxList(propsDict_.lookup("throwingBoxes"));
// 
//     boxes_.setSize(boxList.size());
//     
//     forAll(boxList, b)
//     {
//         const entry& boxI = boxList[b];
//         const dictionary& dict = boxI.dict();
// 
//         vector startPoint = dict.lookup("startPoint");
//         vector endPoint = dict.lookup("endPoint");
//         boxes_[b].resetBoundedBox(startPoint, endPoint);
//     }
//     
//     if (boxList.size() != borderList_.size()) 
//     {
//         FatalError
//             << "Something went wrong with the jamaratControl - setBoundBoxes " << endl
//             << nl << abort(FatalError);  
//         
//     }     
// }


void jamaratControl::controlAfterVelocityII()
{}

void jamaratControl::calculateProperties()
{}

void jamaratControl::output
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

scalar jamaratControl::gaussianDistribution
(
    const scalar& mean,
    const scalar& range
 
)
{
    return mean + 0.25*range*cloud_.rndGen().GaussNormal();
}


} // End namespace Foam

// ************************************************************************* //
