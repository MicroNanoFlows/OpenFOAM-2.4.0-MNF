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

#include "mdDsmcCoupling.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(mdDsmcCoupling, 0);

addToRunTimeSelectionTable(polyCouplingController, mdDsmcCoupling, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
mdDsmcCoupling::mdDsmcCoupling
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict,
    couplingInterface2d &twoDInterfaces,
    couplingInterface3d &threeDInterfaces
)
:
    polyCouplingController(t, molCloud, dict, twoDInterfaces, threeDInterfaces),
    propsDict_(dict.subDict(typeName + "Properties")),
    propsDictSend_(dict.subDict(typeName + "Sending")),
    propsDictRecv_(dict.subDict(typeName + "Receiving")),
    molIds_(),
    output_(false),
#ifdef USE_MUI
    cellCentres_(),
    sendInterfaces_(),
    recvInterfaces_(),
#endif
    sending_(false),
    receiving_(false),
    idList(molCloud_.cP().molIds()),
    rU_(molCloud_.redUnits()),
    fixedBounds_(false),
    fixedRegion_(false),
    fixedRegionMin_(vector::zero),
    fixedRegionMax_(vector::zero),
    fixedBoundMin_(vector::zero),
    fixedBoundMax_(vector::zero),
    fixedBoundNorm_(vector::zero),
    fixedBoundZeroThick_(vector(-1, -1, -1)),
    overlapEnergyLimit_(molCloud_.pot().potentialEnergyLimit()),
    overlapIterations_(100),
    prevMolCount_(0),
    currIteration_(0),
    boundCorr_(0)
{
#ifdef USE_MUI
    //- Determine sending interfaces if defined
    if(propsDictSend_.found("sendingInterfaces"))
    {
        const List<word> interfaces(propsDictSend_.lookup("sendingInterfaces"));

        forAll(interfaces, i)
        {
            //- Find MUI interfaces
            for(size_t j=0; j<threeDInterfaces.interfaces->size(); j++)
            {
                //- If the MUI interface is found then create a copy of its pointer address and store in sendInterfaces_
                if(threeDInterfaces.interfaces->getInterfaceName(j).compare(interfaces[i]) == 0)
                {
                    sendInterfaces_.append(threeDInterfaces.interfaces->getInterface(j));
                    sendInterfaceNames_.append(interfaces[i]); //- Store the receiving interface name
                    break;
                }
            }
        }

        //- Check all interfaces were found
        forAll(sendInterfaces_, i)
        {
        	std::cout << "mdDsmcCoupling::mdDsmcCoupling(): Found 3D MUI coupling interface ("
				      << interfaces[i] << ") to send for domain " << threeDInterfaces.domainName << std::endl;
        }
    }

    //- Determine receiving interfaces if defined
    if(propsDictRecv_.found("receivingInterfaces"))
    {
        const List<word> interfaces(propsDictRecv_.lookup("receivingInterfaces"));

        forAll(interfaces, i)
        {
            //- Find MUI interfaces
            for(size_t j=0; j<threeDInterfaces.interfaces->size(); ++j)
            {
                //- If the MUI interface is found then create a copy of its pointer address and store in sendInterfaces_
                if(threeDInterfaces.interfaces->getInterfaceName(j).compare(interfaces[i]) == 0)
                {
                    recvInterfaces_.append(threeDInterfaces.interfaces->getInterface(j));
                    recvInterfaceNames_.append(interfaces[i]); //- Store the receiving interface name
                    break;
                }
            }
        }

        //- Check all interfaces were found
        forAll(recvInterfaces_, i)
        {
        	std::cout << "mdDsmcCoupling::mdDsmcCoupling(): Found 3D MUI coupling interface ("
				      << interfaces[i] << ") to receive for domain " << threeDInterfaces.domainName << std::endl;
        }

        molId_.setSize(recvInterfaces_.size());
        molHistory_.setSize(recvInterfaces_.size());
    }

    if(sendInterfaces_.size() != 0)
    {
        sending_ = true;
    }

    if(recvInterfaces_.size() != 0)
    {
        receiving_ = true;
    }
#else
    FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
                 << "MUI library not enabled at compilation" << exit(FatalError);
#endif

    writeInTimeDir_ = true;
    writeInCase_ = true;

	selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();
    molNames_ = ids.molIdNames();

    bool regionMinFound = false;
    bool regionMaxFound = false;

    if (propsDict_.found("fixedRegionMin"))
	{
        fixedRegionMin_ = propsDict_.lookup("fixedRegionMin");
        regionMinFound = true;
	}

    if (propsDict_.found("fixedRegionMax"))
    {
        fixedRegionMax_ = propsDict_.lookup("fixedRegionMax");
        regionMaxFound = true;
    }

    if((regionMinFound && !regionMaxFound) || (regionMaxFound && !regionMinFound))
    {
        FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
                      << "Cannot find both fixedRegionMin and fixedRegionMax"
                      << exit(FatalError);
    }
    else
    {
        fixedRegion_ = true;
    }

    bool boundMinFound = false;
    bool boundMaxFound = false;
    bool boundNormFound = false;

    if (propsDict_.found("fixedBoundMin"))
    {
        fixedBoundMin_ = propsDict_.lookup("fixedBoundMin");
        boundMinFound = true;
    }

    if (propsDict_.found("fixedBoundMax"))
    {
        fixedBoundMax_ = propsDict_.lookup("fixedBoundMax");
        boundMaxFound = true;
    }

    if (propsDict_.found("fixedBoundNorm"))
    {
        fixedBoundNorm_ = propsDict_.lookup("fixedBoundNorm");

        if(fixedBoundMax_[0] - fixedBoundMin_[0] == 0.0)
        {
            fixedBoundZeroThick_[0] = 1;
        }

        if(fixedBoundMax_[1] - fixedBoundMin_[1] == 0.0)
        {
            fixedBoundZeroThick_[1] = 1;
        }

        if(fixedBoundMax_[2] - fixedBoundMin_[2] == 0.0)
        {
            fixedBoundZeroThick_[2] = 1;
        }

        if(fixedBoundZeroThick_[0] == -1 && fixedBoundZeroThick_[1] == -1 && fixedBoundZeroThick_[2] == -1)
        {
            FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
                         << "A fixed boundary must have a zero thickness in at least one direction"
                         << exit(FatalError);
        }

        boundNormFound = true;
    }

    if((boundMinFound && !boundMaxFound && !boundNormFound) ||
       (boundMaxFound && !boundMinFound && !boundNormFound) ||
       (boundNormFound && !boundMinFound && !boundMaxFound))
    {
        FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
                     << "Cannot find fixedBoundMin, fixedBoundMax and fixedBoundNorm together"
                     << exit(FatalError);
    }
    else
    {
        fixedBounds_ = true;

        const cellList& cells = mesh_.cells();
        const List<point>& pts = mesh_.points();

        vector meshExtents = mesh_.bounds().max() - mesh_.bounds().min();

        vector boundCorr = meshExtents * 1e-4;

        //- Ensure boundary correction value not larger than 1e-8
        if(boundCorr[0] > 1e-8)
        {
            boundCorr[0] = 1e-8;
        }

        if(boundCorr[1] > 1e-8)
        {
            boundCorr[1] = 1e-8;
        }

        if(boundCorr[2] > 1e-8)
        {
            boundCorr[2] = 1e-8;
        }

        // Pick largest correction value as global
        if(boundCorr[0] > boundCorr[1] && boundCorr[0] > boundCorr[2])
        {
            boundCorr_ = boundCorr[0];
        }

        if(boundCorr[1] > boundCorr[0] && boundCorr[1] > boundCorr[2])
        {
            boundCorr_ = boundCorr[1];
        }

        if(boundCorr[2] > boundCorr[0] && boundCorr[2] > boundCorr[1])
        {
            boundCorr_ = boundCorr[2];
        }

        if(boundCorr[0] == boundCorr[1] && boundCorr[0] == boundCorr[2])
        {
            boundCorr_ = boundCorr[0];
        }

        point cellMin;
        point cellMax;

        //- Determine which cells the fixed boundary intersects
        forAll(cells, cell)
        {
            const labelList& pointList = mesh_.cellPoints(cell);

            cellMin[0] = VGREAT;
            cellMin[1] = VGREAT;
            cellMin[2] = VGREAT;
            cellMax[0] = -VSMALL;
            cellMax[1] = -VSMALL;
            cellMax[2] = -VSMALL;

            forAll(pointList, cellPoint)
            {
                if(pts[pointList[cellPoint]][0] < cellMin[0])
                {
                    cellMin[0] = pts[pointList[cellPoint]][0];
                }

                if(pts[pointList[cellPoint]][0] > cellMax[0])
                {
                    cellMax[0] = pts[pointList[cellPoint]][0];
                }

                if(pts[pointList[cellPoint]][1] < cellMin[1])
                {
                    cellMin[1] = pts[pointList[cellPoint]][1];
                }

                if(pts[pointList[cellPoint]][1] > cellMax[1])
                {
                    cellMax[1] = pts[pointList[cellPoint]][1];
                }

                if(pts[pointList[cellPoint]][2] < cellMin[2])
                {
                    cellMin[2] = pts[pointList[cellPoint]][2];
                }

                if(pts[pointList[cellPoint]][2] > cellMax[2])
                {
                    cellMax[2] = pts[pointList[cellPoint]][2];
                }
            }

            if(fixedBoundZeroThick_[0] == 1) //- fixedBoundMin_ and fixedBoundMax_ are the same in the x
            {
                if(fixedBoundNorm_[0] != 0)
                {
                    scalar boundaryExtend = meshExtents[0] * 0.01;
                    bool test = false;

                    if((fixedBoundMin_[0]-boundaryExtend) >= cellMin[0] && (fixedBoundMin_[0]-boundaryExtend) <= cellMax[0])
                    {
                        test = true;
                    }

                    if((fixedBoundMin_[0]+boundaryExtend) >= cellMin[0] && (fixedBoundMin_[0]+boundaryExtend) <= cellMax[0])
                    {
                        test = true;
                    }

                    if(test)
                    {
                        if((cellMin[1] >= fixedBoundMin_[1] && cellMax[1] <= fixedBoundMax_[1]) &&
                           (cellMin[2] >= fixedBoundMin_[2] && cellMax[2] <= fixedBoundMax_[2]))
                        {
                            if(findIndex(intersectingCells_, cell) == -1)
                            {
                                intersectingCells_.append(cell);
                            }
                        }
                    }
                }
                else
                {
                    FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
                                 << "Fixed boundary zero thickness in x direction but normal value zero"
                                 << exit(FatalError);
                }
            }

            if(fixedBoundZeroThick_[1] == 1) //- fixedBoundMin_ and fixedBoundMax_ are the same in the y
            {
               if(fixedBoundNorm_[1] != 0)
               {
                   scalar boundaryExtend = meshExtents[1] * 0.01;
                   bool test = false;

                   if((fixedBoundMin_[1]-boundaryExtend) >= cellMin[1] && (fixedBoundMin_[1]-boundaryExtend) <= cellMax[1])
                   {
                       test = true;
                   }

                   if((fixedBoundMin_[1]+boundaryExtend) >= cellMin[1] && (fixedBoundMin_[1]+boundaryExtend) <= cellMax[1])
                   {
                       test = true;
                   }

                   if(test)
                   {
                       if((cellMin[0] >= fixedBoundMin_[0] && cellMax[0] <= fixedBoundMax_[0]) &&
                          (cellMin[2] >= fixedBoundMin_[2] && cellMax[2] <= fixedBoundMax_[2]))
                       {
                           if(findIndex(intersectingCells_, cell) == -1)
                           {
                               intersectingCells_.append(cell);
                           }
                       }
                   }
               }
               else
               {
                   FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
                                << "Fixed boundary zero thickness in y direction but normal value zero"
                                << exit(FatalError);
               }
            }

            if(fixedBoundZeroThick_[2] == 1) //- fixedBoundMin_ and fixedBoundMax_ are the same in the z
            {
               if(fixedBoundNorm_[2] != 0)
               {
                   scalar boundaryExtend = meshExtents[2] * 0.01;
                   bool test = false;

                   if((fixedBoundMin_[2]-boundaryExtend) >= cellMin[2] && (fixedBoundMin_[2]-boundaryExtend) <= cellMax[2])
                   {
                       test = true;
                   }

                   if((fixedBoundMin_[2]+boundaryExtend) >= cellMin[2] && (fixedBoundMin_[2]+boundaryExtend) <= cellMax[2])
                   {
                       test = true;
                   }

                   if(test)
                   {
                       if((cellMin[0] >= fixedBoundMin_[0] && cellMax[0] <= fixedBoundMax_[0]) &&
                          (cellMin[1] >= fixedBoundMin_[1] && cellMax[1] <= fixedBoundMax_[1]))
                       {
                           if(findIndex(intersectingCells_, cell) == -1)
                           {
                               intersectingCells_.append(cell);
                           }
                       }
                   }
               }
               else
               {
                   FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
                                << "Fixed boundary zero thickness in z direction but normal value zero"
                                << exit(FatalError);
               }
            }
        }

        if(intersectingCells_.size() > 0)
        {
            std::cout << "Fixed boundary intersecting cell count: " << intersectingCells_.size() << std::endl;
        }
        else
        {
            //- Determine if fixed boundary falls within mesh bounds
            point meshMin(VGREAT, VGREAT, VGREAT);
            point meshMax(-VSMALL, -VSMALL, -VSMALL);

            const pointField& meshPoints = mesh_.points();

            forAll(meshPoints, pts)
            {
                if(meshPoints[pts][0] < meshMin[0])
                {
                    meshMin[0] = meshPoints[pts][0];
                }

                if(meshPoints[pts][1] < meshMin[1])
                {
                    meshMin[1] = meshPoints[pts][1];
                }

                if(meshPoints[pts][2] < meshMin[2])
                {
                    meshMin[2] = meshPoints[pts][2];
                }

                if(meshPoints[pts][0] > meshMax[0])
                {
                    meshMax[0] = meshPoints[pts][0];
                }

                if(meshPoints[pts][1] > meshMax[1])
                {
                    meshMax[1] = meshPoints[pts][1];
                }

                if(meshPoints[pts][2] > meshMax[2])
                {
                    meshMax[2] = meshPoints[pts][2];
                }
            }

            vector meshHalfWidth(((meshMax[0] - meshMin[0]) * 0.5),
                                 ((meshMax[1] - meshMin[1]) * 0.5),
                                 ((meshMax[2] - meshMin[2]) * 0.5));
            vector fixedBoundHalfWidth(((fixedBoundMax_[0] - fixedBoundMin_[0]) * 0.5),
                                       ((fixedBoundMax_[1] - fixedBoundMin_[1]) * 0.5),
                                       ((fixedBoundMax_[2] - fixedBoundMin_[2]) * 0.5));
            point meshCentre(meshMin[0] + meshHalfWidth[0],
                             meshMin[1] + meshHalfWidth[1],
                             meshMin[2] + meshHalfWidth[2]);
            point fixedBoundCentre(fixedBoundMin_[0] + fixedBoundHalfWidth[0],
                                   fixedBoundMin_[1] + fixedBoundHalfWidth[1],
                                   fixedBoundMin_[2] + fixedBoundHalfWidth[2]);

            bool noOverlap = false;

            if ((std::fabs(meshCentre[0] - fixedBoundCentre[0]) > (meshHalfWidth[0] + fixedBoundHalfWidth[0])) ||
               (std::fabs(meshCentre[1] - fixedBoundCentre[1]) > (meshHalfWidth[1] + fixedBoundHalfWidth[1])) ||
               (std::fabs(meshCentre[2] - fixedBoundCentre[2]) > (meshHalfWidth[2] + fixedBoundHalfWidth[2])))
            {
                noOverlap = true;
            }

            //- There is an overlap between the fixed boundary and the local mesh so should have found at least 1 intersecting cell
            if(!noOverlap)
            {
                FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
                             << "Fixed boundary defined but no intersecting cells found"
                             << exit(FatalError);
            }
        }
    }

    if (propsDict_.found("potentialEnergyInsertLimit"))
    {
        overlapEnergyLimit_ = readScalar(propsDict_.lookup("potentialEnergyInsertLimit"));
    }

    if (propsDict_.found("insertIterations"))
    {
        overlapIterations_ = readLabel(propsDict_.lookup("insertIterations"));
    }

    if (propsDict_.found("output"))
    {
        output_ = Switch(propsDict_.lookup("output"));
    }

    if(sending_ || receiving_)
	{
		refLength_ = threeDInterfaces.refLength; //- Store the reference length
		refTime_ = threeDInterfaces.refTime; //- Store the reference time

		oneOverRefLength_ = 1.0 / refLength_;
		oneOverRefTime_ = 1.0 / refTime_;
#ifdef USE_MUI
		//Initialise exact time sampler for MUI
		chrono_sampler = new mui::chrono_sampler_exact3d();
#endif
	}

    refVelocity_ = rU_.refVelocity();
    oneOverRefVelocity_ = 1.0 / refVelocity_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mdDsmcCoupling::~mdDsmcCoupling()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mdDsmcCoupling::initialConfiguration()
{
	if(receiving_)
    {
		receiveCoupledRegion(true); // Receive ghost molecules in coupled regions at time = startTime and commit time=1 to release other side
    }
}

void mdDsmcCoupling::controlAfterMove(label stage)
{
	if(stage == 1)
	{
        currIteration_++; //- Increment the current iteration

        if(sending_)
        {
            findCoupledMolecules(); //- Find, collate and delete any molecules that have passed a coupling boundary
            sendCoupledMolecules(); // Send any molecules deleted by coupling boundary (non-blocking)
        }
	}
	else if (stage == 2)
	{
	    if(receiving_)
        {
            receiveCoupledParcels(); // Receive any molecules coupling boundary (blocking)
        }
	}
	else if (stage == 3)
	{
		if(receiving_)
        {
            receiveCoupledRegion(false); // Receive ghost molecules in coupled region (blocking)
        }
    }
	else if (stage == 4)
	{
	    if(receiving_)
	    {
	        insertCoupledMolecules(); // Insert any coupling boundary molecules from stage 2
	    }
	}
}

void mdDsmcCoupling::controlAfterForces()
{
    if(sending_)
    {
        sendCoupledRegionForces(); // Send the forces acting on molecules in the coupling region(s) (non-blocking)
    }
}

void mdDsmcCoupling::receiveCoupledRegion(bool init)
{
#ifdef USE_MUI
	label molCount = 0;
	List<std::vector<mui::point3d> > rcvPoints(recvInterfaces_.size());
	List<std::vector<std::string> > rcvMolType(recvInterfaces_.size());
	List<std::vector<label> > rcvMolId(recvInterfaces_.size());
    List<std::vector<scalar> > rcvVelX(recvInterfaces_.size());
    List<std::vector<scalar> > rcvVelY(recvInterfaces_.size());
    List<std::vector<scalar> > rcvVelZ(recvInterfaces_.size());

	// Iterate through all receiving interfaces for this controller and extract a points list
	forAll(recvInterfaces_, iface)
	{
        //- Extract a list of all molecule locations received from other solver through this interface
        rcvPoints[iface] = recvInterfaces_[iface]->fetch_points<std::string>("type_region", currIteration_, *chrono_sampler);

        if(rcvPoints[iface].size() > 0)
        {
            //- Extract a list of all molecule change status values received from other solver through this interface
            rcvMolType[iface] = recvInterfaces_[iface]->fetch_values<std::string>("type_region", currIteration_, *chrono_sampler);

            //- Extract a list of all molecule Id's received from other solver through this interface
            rcvMolId[iface] = recvInterfaces_[iface]->fetch_values<label>("id_region", currIteration_, *chrono_sampler);

            //- Extract a list of all molecule velocities received from other solver through this interface
            rcvVelX[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_x_region", currIteration_, *chrono_sampler);
            rcvVelY[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_y_region", currIteration_, *chrono_sampler);
            rcvVelZ[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_z_region", currIteration_, *chrono_sampler);
        }
	}

	//- Go through received values and find any that are not of the type set to be received
    forAll(rcvMolType, ifacepts)
    {
        if(rcvMolType[ifacepts].size() > 0)
        {
            std::vector<std::string>::iterator rcvMolTypeIt;
            std::vector<mui::point3d>::iterator rcvPointsIt = rcvPoints[ifacepts].begin();
            std::vector<label>::iterator rcvMolIdIt = rcvMolId[ifacepts].begin();
            std::vector<scalar>::iterator rcvVelXIt = rcvVelX[ifacepts].begin();
            std::vector<scalar>::iterator rcvVelYIt = rcvVelY[ifacepts].begin();
            std::vector<scalar>::iterator rcvVelZIt = rcvVelZ[ifacepts].begin();

            for (rcvMolTypeIt = rcvMolType[ifacepts].begin(); rcvMolTypeIt != rcvMolType[ifacepts].end(); rcvMolTypeIt++) {
                const label molId = findIndex(molNames_, *rcvMolTypeIt);

                if(molId == -1) //- molId not found in local list as one to receive so store it as one to remove from lists
                {
                    rcvMolType[ifacepts].erase(rcvMolTypeIt--);
                    rcvPoints[ifacepts].erase(rcvPointsIt--);
                    rcvMolId[ifacepts].erase(rcvMolIdIt--);
                    rcvVelX[ifacepts].erase(rcvVelXIt--);
                    rcvVelY[ifacepts].erase(rcvVelYIt--);
                    rcvVelZ[ifacepts].erase(rcvVelZIt--);
                }

                rcvPointsIt++;
                rcvMolIdIt++;
                rcvVelXIt++;
                rcvVelYIt++;
                rcvVelZIt++;
            }
        }
    }

	//- Insert/update the ghost molecules
	forAll(rcvPoints, ifacepts)
	{
		if(rcvPoints[ifacepts].size() > 0)
		{
		    bool newList = false;
			
			//If list size has changed then treat this as a new list
			if(prevMolCount_ != rcvMolId[ifacepts].size())
			{
				if(molId_[ifacepts].size() > 0)
				{
					forAll(molHistory_[ifacepts], mol)
					{
					    molCloud_.removeMolFromCellOccupancy(molHistory_[ifacepts][mol]->origId(), molHistory_[ifacepts][mol]->cell());
					    molCloud_.deleteParticle(*molHistory_[ifacepts][mol]);
					}

					molId_[ifacepts].clear();
					molHistory_[ifacepts].clear();
				}

				molId_[ifacepts].setSize(rcvMolId[ifacepts].size(), -1);
				molHistory_[ifacepts].setSize(rcvMolId[ifacepts].size(), NULL);
				newList = true;
			}

            for (size_t pts = 0; pts < rcvPoints[ifacepts].size(); pts++)
            {
                const label molId = findIndex(molNames_, rcvMolType[ifacepts][pts]);

                vector velocity;
                velocity[0] = rcvVelX[ifacepts][pts] * oneOverRefVelocity_;
                velocity[1] = rcvVelY[ifacepts][pts] * oneOverRefVelocity_;
                velocity[2] = rcvVelZ[ifacepts][pts] * oneOverRefVelocity_;

                point checkedPosition((rcvPoints[ifacepts][pts][0] * refLength_), (rcvPoints[ifacepts][pts][1] * refLength_), (rcvPoints[ifacepts][pts][2] * refLength_));

                if(fixedRegion_)
                {
                    if(checkedPosition[0] <= fixedRegionMin_[0])
                    {
                        checkedPosition[0] = fixedRegionMin_[0] + boundCorr_;
                    }

                    if(checkedPosition[0] >= fixedRegionMax_[0])
                    {
                        checkedPosition[0] = fixedRegionMax_[0] - boundCorr_;
                    }

                    if(checkedPosition[1] <= fixedRegionMin_[1])
                    {
                        checkedPosition[1] = fixedRegionMin_[1] + boundCorr_;
                    }

                    if(checkedPosition[1] >= fixedRegionMax_[1])
                    {
                        checkedPosition[1] = fixedRegionMax_[1] - boundCorr_;
                    }

                    if(checkedPosition[2] <= fixedRegionMin_[2])
                    {
                        checkedPosition[2] = fixedRegionMin_[2] + boundCorr_;
                    }

                    if(checkedPosition[2] >= fixedRegionMax_[2])
                    {
                        checkedPosition[2] = fixedRegionMax_[2] - boundCorr_;
                    }
                }

                if(newList) //This is a completely new list so all molecules to be inserted
                {
                    molId_[ifacepts][pts] = rcvMolId[ifacepts][pts];
                    molHistory_[ifacepts][pts] = insertMolecule(checkedPosition, molIds_[molId], true, velocity);
                    molCount++;
                }
                else //This is not a completely new list so just check for individual molecule changes
                {
                    if(molId_[ifacepts][pts] != rcvMolId[ifacepts][pts]) //This molecule has changed in the list since the last time
                    {
                        molId_[ifacepts][pts] = rcvMolId[ifacepts][pts]; // Update ID history

                        molCloud_.removeMolFromCellOccupancy(molHistory_[ifacepts][pts]->origId(), molHistory_[ifacepts][pts]->cell());
                        molCloud_.deleteParticle(*molHistory_[ifacepts][pts]); //Delete the molecule in this list position

                        molHistory_[ifacepts][pts] = insertMolecule(checkedPosition, molIds_[molId], true, velocity); //Insert the new molecule

                        molCount++;
                    }
                    else //This molecule is the same as received last time, so just update it's properties
                    {
                        molHistory_[ifacepts][pts]->position()[0] = checkedPosition[0];
                        molHistory_[ifacepts][pts]->position()[1] = checkedPosition[1];
                        molHistory_[ifacepts][pts]->position()[2] = checkedPosition[2];

                        molHistory_[ifacepts][pts]->v()[0] = velocity[0];
                        molHistory_[ifacepts][pts]->v()[1] = velocity[1];
                        molHistory_[ifacepts][pts]->v()[2] = velocity[2];

                        molCount++;
                    }
                }
            }
		}
	}

	if(init)
	{
		forAll(recvInterfaces_, iface)
		{
		    recvInterfaces_[iface]->commit(currIteration_);
		}
	}

	if(molCount > 0)
	{
        if(molCount != prevMolCount_)
        {
            if (Pstream::parRun())
            {
                if(Pstream::master())
                {
                    if(prevMolCount_ > 0)
                    {
                        std::cout << "Number of molecules in coupled region now equals " << molCount << " (previously " << prevMolCount_ << ")" << std::endl;
                    }
                    else
                    {
                        std::cout << "Number of molecules in coupled region now equals " << molCount << std::endl;
                    }
                }
                else
                {
                    if(prevMolCount_ > 0)
                    {
                        std::cout << "[" << time_.time().value() << "s] Number of molecules in coupled region now equals " << molCount << " (previously " << prevMolCount_ << ")" << std::endl;
                    }
                    else
                    {
                        std::cout << "[" << time_.time().value() << "s] Number of molecules in coupled region now equals " << molCount << std::endl;
                    }
                }
            }
            else
            {
                if(prevMolCount_ > 0)
                {
                    std::cout << "Number of molecules in coupled region now equals " << molCount << " (previously " << prevMolCount_ << ")" << std::endl;
                }
                else
                {
                    std::cout << "Number of molecules in coupled region now equals " << molCount << std::endl;
                }
            }

            prevMolCount_ = molCount;
        }
	}
#endif
}

void mdDsmcCoupling::sendCoupledRegionForces()
{
#ifdef USE_MUI
    polyMolecule* molecule = NULL;

    // Iterate through all sending interfaces for this controller
    forAll(sendInterfaces_, iface)
    {
        forAll(molHistory_[iface], mol) // Iterate through molecules received by this controller
        {
            molecule = molHistory_[iface][mol];

            //- Determine whether parcel is of a type set to send by this controller
            const label typeIndex = findIndex(molNames_, molCloud_.cP().molIds()[molecule->id()]);

            if(typeIndex != -1)
            {
                // Get the molecule centre
                mui::point3d molCentre;
                molCentre[0] = molecule->position()[0] * oneOverRefLength_;
                molCentre[1] = molecule->position()[1] * oneOverRefLength_;
                molCentre[2] = molecule->position()[2] * oneOverRefLength_;

                // Push molecule type
                sendInterfaces_[iface]->push("type_region", molCentre, static_cast<std::string>(molNames_[typeIndex]));

                // Push molecule ID from receive history
                sendInterfaces_[iface]->push("id_region", molCentre, static_cast<label>(molId_[iface][mol]));

                vector siteForcesAccum(vector::zero);

                forAll(molecule->siteForces(), s)
                {
                    siteForcesAccum[0] += molecule->siteForces()[s][0];
                    siteForcesAccum[1] += molecule->siteForces()[s][1];
                    siteForcesAccum[2] += molecule->siteForces()[s][2];
                }

                // Push the molecule site forces to the interface
                sendInterfaces_[iface]->push("force_x_region", molCentre, siteForcesAccum[0] * rU_.refForce());
                sendInterfaces_[iface]->push("force_y_region", molCentre, siteForcesAccum[1] * rU_.refForce());
                sendInterfaces_[iface]->push("force_z_region", molCentre, siteForcesAccum[2] * rU_.refForce());
            }
        }
    }

    forAll(sendInterfaces_, iface)
    {
        // Commit (transmit) values to the MUI interface
        sendInterfaces_[iface]->commit(currIteration_);
    }
#endif
}

void mdDsmcCoupling::findCoupledMolecules()
{
    DynamicList<polyMolecule*> molsToRemove;

    forAll(intersectingCells_, cell)
    {
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[intersectingCells_[cell]];

        forAll(molsInCell, molecule)
        {
            const word& molType = molCloud_.cP().molIds()[molsInCell[molecule]->id()];
            const label typeIndex = findIndex(molNames_, molType);

            //- Only delete and store particles that are of the coupled type
            if(typeIndex != -1)
            {
                if(!molsInCell[molecule]->ghost() && !molsInCell[molecule]->frozen())
                {
                    bool removeMolecule = false;

                    if(fixedBoundZeroThick_[0] == 1)
                    {
                        if(fixedBoundNorm_[0] > 0) //- Boundary is positive facing in the x
                        {
                            if(molsInCell[molecule]->position()[0] <= fixedBoundMin_[0])
                            {
                                removeMolecule = true;
                            }
                        }
                        else if (fixedBoundNorm_[0] < 0) //- Boundary is negative facing in the x
                        {
                            if(molsInCell[molecule]->position()[0] >= fixedBoundMax_[0])
                            {
                                removeMolecule = true;
                            }
                        }
                    }

                    if(fixedBoundZeroThick_[1] == 1)
                    {
                        if(fixedBoundNorm_[1] > 0) //- Boundary is positive facing in the y
                        {
                            if(molsInCell[molecule]->position()[1] <= fixedBoundMin_[1])
                            {
                                removeMolecule = true;
                            }
                        }
                        else if (fixedBoundNorm_[1] < 0) //- Boundary is negative facing in the x
                        {
                            if(molsInCell[molecule]->position()[1] >= fixedBoundMax_[1])
                            {
                                removeMolecule = true;
                            }
                        }
                    }

                    if(fixedBoundZeroThick_[2] == 1)
                    {
                        if(fixedBoundNorm_[2] > 0) //- Boundary is positive facing in the z
                        {
                            if(molsInCell[molecule]->position()[2] <= fixedBoundMin_[2])
                            {
                                removeMolecule = true;
                            }
                        }
                        else if (fixedBoundNorm_[2] < 0) //- Boundary is negative facing in the x
                        {
                            if(molsInCell[molecule]->position()[2] >= fixedBoundMax_[2])
                            {
                                removeMolecule = true;
                            }
                        }
                    }

                    if(removeMolecule)
                    {
                        //- Store the required details of the molecule
                        coupledMolecule newMolToSend;
                        newMolToSend.molType = molCloud_.cP().molIds()[molsInCell[molecule]->id()];
                        newMolToSend.position = molsInCell[molecule]->position();

                        if(fixedBounds_)
                        {
                            if(fixedBoundZeroThick_[0] == 1) //- Boundary has zero thickness in the x
                            {
                                newMolToSend.position[0] = fixedBoundMin_[0];
                            }
                            else
                            {
                                if(newMolToSend.position[0] <= fixedBoundMin_[0])
                                {
                                    newMolToSend.position[0] = fixedBoundMin_[0];
                                }

                                if(newMolToSend.position[0] >= fixedBoundMax_[0])
                                {
                                    newMolToSend.position[0] = fixedBoundMax_[0];
                                }
                            }

                            if(fixedBoundZeroThick_[1] == 1) //- Boundary has zero thickness in the y
                            {
                                newMolToSend.position[1] = fixedBoundMin_[1];
                            }
                            else
                            {
                                if(newMolToSend.position[1] <= fixedBoundMin_[1])
                                {
                                    newMolToSend.position[1] = fixedBoundMin_[1];
                                }

                                if(newMolToSend.position[1] >= fixedBoundMax_[1])
                                {
                                    newMolToSend.position[1] = fixedBoundMax_[1];
                                }
                            }

                            if(fixedBoundZeroThick_[2] == 1) //- Boundary has zero thickness in the z
                            {
                                newMolToSend.position[2] = fixedBoundMin_[2];
                            }
                            else
                            {
                                if(newMolToSend.position[2] <= fixedBoundMin_[2])
                                {
                                    newMolToSend.position[2] = fixedBoundMin_[2];
                                }

                                if(newMolToSend.position[2] >= fixedBoundMax_[2])
                                {
                                    newMolToSend.position[2] = fixedBoundMax_[2];
                                }
                            }
                        }

                        newMolToSend.velocity = molsInCell[molecule]->v();
                        molsToSend_.append(newMolToSend);

                        //- Store that this molecule needs to be removed
                        molsToRemove.append(molsInCell[molecule]);
                    }
                }
            }
        }
    }

    //- Iterate through all parcels to be removed
    forAll(molsToRemove, molecule)
    {
        //- Delete molecule from cellOccupancy (before deleting it from cloud)
        molCloud_.removeMolFromCellOccupancy(molsToRemove[molecule]->origId(), molsToRemove[molecule]->cell());
        molCloud_.deleteParticle(*molsToRemove[molecule]);
    }
}

void mdDsmcCoupling::sendCoupledMolecules()
{
#ifdef USE_MUI
    if(molsToSend_.size() != 0)
	{
		forAll(molsToSend_, mols)
		{
            forAll(sendInterfaces_, iface)
            {
                // Get the molecule centre
                mui::point3d molCentre;
                molCentre[0] = molsToSend_[mols].position[0] * oneOverRefLength_;
                molCentre[1] = molsToSend_[mols].position[1] * oneOverRefLength_;
                molCentre[2] = molsToSend_[mols].position[2] * oneOverRefLength_;

                // Push molecule type
                sendInterfaces_[iface]->push("type_bound", molCentre, static_cast<std::string>(molsToSend_[mols].molType));

                // Push molecule velocity
                sendInterfaces_[iface]->push("vel_x_bound", molCentre, molsToSend_[mols].velocity[0] * rU_.refVelocity());
                sendInterfaces_[iface]->push("vel_y_bound", molCentre, molsToSend_[mols].velocity[1] * rU_.refVelocity());
                sendInterfaces_[iface]->push("vel_z_bound", molCentre, molsToSend_[mols].velocity[2] * rU_.refVelocity());

                if (Pstream::parRun())
                {
                    if(Pstream::master())
                    {
                        std::cout << "Coupling boundary molecule pushed at: [" << molCentre[0] << "," << molCentre[1] << "," << molCentre[2] << "]" << std::endl;
                    }
                    else
                    {
                        std::cout << "[" << time_.time().value() << "s] Coupling boundary molecule pushed at: [" << molCentre[0] << "," << molCentre[1] << "," << molCentre[2] << "]" << std::endl;
                    }
                }
                else
                {
                    std::cout << "Coupling boundary molecule pushed at: [" << molCentre[0] << "," << molCentre[1] << "," << molCentre[2] << "]" << std::endl;
                }
            }
		}

		//- Clear the sent molecules
		molsToSend_.clear();
	}

    forAll(sendInterfaces_, iface)
    {
    	// Commit values to the coupling interface
	  	sendInterfaces_[iface]->commit(currIteration_);
    }
#endif
}

void mdDsmcCoupling::receiveCoupledParcels()
{
#ifdef USE_MUI
	List<std::vector<mui::point3d> > rcvPoints(recvInterfaces_.size());
	List<std::vector<std::string> > rcvMolType(recvInterfaces_.size());
    List<std::vector<scalar> > rcvVelX(recvInterfaces_.size());
    List<std::vector<scalar> > rcvVelY(recvInterfaces_.size());
    List<std::vector<scalar> > rcvVelZ(recvInterfaces_.size());

	// Iterate through all receiving interfaces for this controller and extract a points list for each molecule type handled
	forAll(recvInterfaces_, iface)
	{
        //- Extract a list of all molecule locations
        rcvPoints[iface] = recvInterfaces_[iface]->fetch_points<std::string>("type_bound", currIteration_, *chrono_sampler);

        if(rcvPoints[iface].size() > 0)
        {
            //- Extract a list of all molecule types
            rcvMolType[iface] = recvInterfaces_[iface]->fetch_values<std::string>("type_bound", currIteration_, *chrono_sampler);

            //- Extract a list of all molecule velocities received from other solver through this interface
            rcvVelX[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_x_bound", currIteration_, *chrono_sampler);
            rcvVelY[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_y_bound", currIteration_, *chrono_sampler);
            rcvVelZ[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_z_bound", currIteration_, *chrono_sampler);
        }
	}

	//- Go through received values and find any that are not of the type set to be received
    forAll(rcvMolType, ifacepts)
    {
        if(rcvMolType[ifacepts].size() > 0)
        {
            std::vector<std::string>::iterator rcvMolTypeIt;
            std::vector<mui::point3d>::iterator rcvPointsIt = rcvPoints[ifacepts].begin();
            std::vector<scalar>::iterator rcvVelXIt = rcvVelX[ifacepts].begin();
            std::vector<scalar>::iterator rcvVelYIt = rcvVelY[ifacepts].begin();
            std::vector<scalar>::iterator rcvVelZIt = rcvVelZ[ifacepts].begin();

            for (rcvMolTypeIt = rcvMolType[ifacepts].begin(); rcvMolTypeIt != rcvMolType[ifacepts].end(); rcvMolTypeIt++) {
                const label molId = findIndex(molNames_, *rcvMolTypeIt);

                if(molId == -1) //- molId not found in local list as one to receive so store it as one to remove from lists
                {
                    rcvMolType[ifacepts].erase(rcvMolTypeIt--);
                    rcvPoints[ifacepts].erase(rcvPointsIt--);
                    rcvVelX[ifacepts].erase(rcvVelXIt--);
                    rcvVelY[ifacepts].erase(rcvVelYIt--);
                    rcvVelZ[ifacepts].erase(rcvVelZIt--);
                }

                rcvPointsIt++;
                rcvVelXIt++;
                rcvVelYIt++;
                rcvVelZIt++;
            }
        }
    }

    // Iterate through all receiving interfaces for this controller
	forAll(recvInterfaces_, ifacepts)
	{
        if(rcvPoints[ifacepts].size() > 0)
        {
            for (size_t pts = 0; pts < rcvPoints[ifacepts].size(); pts++)
            {
                vector velocity;
                velocity[0] = rcvVelX[ifacepts][pts] * oneOverRefVelocity_;
                velocity[1] = rcvVelY[ifacepts][pts] * oneOverRefVelocity_;
                velocity[2] = rcvVelZ[ifacepts][pts] * oneOverRefVelocity_;

                point checkedPosition(rcvPoints[ifacepts][pts][0] * refLength_, rcvPoints[ifacepts][pts][1] * refLength_, rcvPoints[ifacepts][pts][2] * refLength_);

                if(fixedBounds_)
                {
                    if(fixedBoundZeroThick_[0] == 1) //- Boundary has zero thickness in the x
                    {
                        checkedPosition[0] = fixedBoundMin_[0] + (fixedBoundNorm_[0] * boundCorr_);
                    }
                    else
                    {
                        if(checkedPosition[0] <= fixedBoundMin_[0])
                        {
                            checkedPosition[0] = fixedBoundMin_[0] + boundCorr_;
                        }

                        if(checkedPosition[0] >= fixedBoundMax_[0])
                        {
                            checkedPosition[0] = fixedBoundMax_[0] - boundCorr_;
                        }
                    }

                    if(fixedBoundZeroThick_[1] == 1) //- Boundary has zero thickness in the y
                    {
                        checkedPosition[1] = fixedBoundMin_[1] + (fixedBoundNorm_[1] * boundCorr_);
                    }
                    else
                    {
                        if(checkedPosition[1] <= fixedBoundMin_[1])
                        {
                            checkedPosition[1] = fixedBoundMin_[1] + boundCorr_;
                        }

                        if(checkedPosition[1] >= fixedBoundMax_[1])
                        {
                            checkedPosition[1] = fixedBoundMax_[1] - boundCorr_;
                        }
                    }

                    if(fixedBoundZeroThick_[2] == 1) //- Boundary has zero thickness in the z
                    {
                        checkedPosition[2] = fixedBoundMin_[2] + (fixedBoundNorm_[2] * boundCorr_);
                    }
                    else
                    {
                        if(checkedPosition[2] <= fixedBoundMin_[2])
                        {
                            checkedPosition[2] = fixedBoundMin_[2] + boundCorr_;
                        }

                        if(checkedPosition[2] >= fixedBoundMax_[2])
                        {
                            checkedPosition[2] = fixedBoundMax_[2] - boundCorr_;
                        }
                    }
                }

                if (Pstream::parRun())
                {
                    if(Pstream::master())
                    {
                        std::cout << "Coupling boundary parcel received at: [" << checkedPosition[0] << "," << checkedPosition[1] << "," << checkedPosition[2] << "]" << std::endl;
                    }
                    else
                    {
                        std::cout << "[" << time_.time().value() << "s] Coupling boundary parcel received at: [" << checkedPosition[0] << "," << checkedPosition[1] << "," << checkedPosition[2] << "]" << std::endl;
                    }
                }
                else
                {
                    std::cout << "Coupling boundary parcel received at: [" << checkedPosition[0] << "," << checkedPosition[1] << "," << checkedPosition[2] << "]" << std::endl;
                }

                coupledMolecule newMol;
                newMol.molType = rcvMolType[ifacepts][pts];
                newMol.position = checkedPosition;
                newMol.velocity = velocity;

                molsReceived_.append(newMol);
            }
        }
	}
#endif
}

bool mdDsmcCoupling::insertCoupledMolecules()
{
    bool molInserted = false;
#ifdef USE_MUI
    if(molsReceived_.size() > 0)
    {
        forAll(molsReceived_, mol)
        {
            const label molId = findIndex(molNames_, molsReceived_[mol].molType);

            point molPosition = molsReceived_[mol].position;

            if(insertMolecule(molPosition, molIds_[molId], false, molsReceived_[mol].velocity) != NULL)
            {
                if (Pstream::parRun())
                {
                    if(Pstream::master())
                    {
                        std::cout << "Coupling boundary parcel inserted at: [" << molPosition[0] << "," << molPosition[1] << "," << molPosition[2] << "]" << std::endl;
                    }
                    else
                    {
                        std::cout << "[" << time_.time().value() << "s] Coupling boundary parcel inserted at: [" << molPosition[0] << "," << molPosition[1] << "," << molPosition[2] << "]" << std::endl;
                    }
                }
                else
                {
                    std::cout << "Coupling boundary parcel inserted at: [" << molPosition[0] << "," << molPosition[1] << "," << molPosition[2] << "]" << std::endl;
                }

                molInserted = true;
            }
        }

        molsReceived_.clear();
    }
#endif
    return molInserted;
}

void mdDsmcCoupling::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    if(output_)
    {
        label singleStepNMols = molCloud_.size();
        scalar singleStepTotalLinearKE = 0.0;

        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for(; mol != molCloud_.end(); ++mol)
        {
            label molId = mol().id();

            scalar molMass(molCloud_.cP().mass(molId));

            singleStepTotalLinearKE += 0.5*molMass*magSqr(mol().v());
        }

        if (Pstream::parRun())
        {
            reduce(singleStepTotalLinearKE, sumOp<scalar>());
            reduce(singleStepNMols, sumOp<label>());
        }

        if(Pstream::master())
        {
            OFstream os1(timePath/"avg_lin_KE");
            os1 << "Time " << time_.time().value() << endl;
            os1 << singleStepTotalLinearKE/singleStepNMols << endl;
        }
    }
}

void mdDsmcCoupling::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateCouplingControllerProperties(newDict);
    propsDict_ = newDict.subDict(typeName + "Properties");
    propsDictSend_ = newDict.subDict(typeName + "Sending");
	propsDictRecv_ = newDict.subDict(typeName + "Receiving");
}

void mdDsmcCoupling::barrier()
{
	forAll(sendInterfaces_, iface)
	{
		sendInterfaces_[iface]->barrier(currIteration_);
	}
}

void mdDsmcCoupling::barrier(label iteration)
{
	forAll(sendInterfaces_, iface)
	{
		sendInterfaces_[iface]->barrier(iteration);
	}
}

void mdDsmcCoupling::barrier(label iteration, label interface)
{
	sendInterfaces_[interface]->barrier(iteration);
}

void mdDsmcCoupling::forget()
{
	forAll(recvInterfaces_, iface)
	{
		recvInterfaces_[iface]->forget(currIteration_, true);
	}
}

void mdDsmcCoupling::forget(label iteration)
{
	forAll(recvInterfaces_, iface)
	{
		recvInterfaces_[iface]->forget(iteration, true);
	}
}

void mdDsmcCoupling::forget(label iteration, label interface)
{
	recvInterfaces_[interface]->forget(iteration, true);
}

polyMolecule* mdDsmcCoupling::insertMolecule
(
    point& position,
    const label& id,
    const bool ghost,
    vector& velocity
)
{
    label cell = mesh_.findCell(position);
    
    if(cell != -1)
    {
        label tetFace = -1;
        label tetPt = -1;

        mesh_.findCellFacePt
        (
            position,
            cell,
            tetFace,
            tetPt
        );

    	point specialPosition(vector::zero);

        label special = 0;

        if (ghost)
        {
            specialPosition = position;
            special = polyMolecule::SPECIAL_GHOST;
        }

        vector pi = vector::zero;

        tensor Q = I;

        polyMolecule* newMol = molCloud_.createMolecule
        (
            position,
            cell,
            tetFace,
            tetPt,
            Q,
            velocity,
            vector::zero,
            pi,
            vector::zero,
            specialPosition,
            special,
            id,
            1.0,
            molCloud_.getTrackingNumber()
        );

        molCloud_.insertMolInCellOccupancy(newMol);

        if(!ghost) // No need to perform overlap check for ghost molecules
        {
            polyMolecule* overlapMol = NULL;
            overlapMol = checkForOverlaps(newMol, overlapEnergyLimit_);
            label iterCount = 0;

            if(overlapMol != NULL) // Check for any overlaps that will exceed pre-determined energy limit for interactions
            {
                std::cout << "mdDsmcCoupling::insertMolecule(): Molecule insertion would create overlap, finding new insertion position" << std::endl;

                vector startPos = position;
                vector overLap(vector::zero);
                scalar overLapMag = 0;
                vector overlapNorm(vector::zero);
                scalar perturbPerc = 0.5;
                scalar perturbDistance = 0;
                bool trunkX = false, trunkY = false, trunkZ = false;
                label stuckIter = 0;
                label randomNorm = -1; // Set at -1 so norm calculated using geometry initially
                label randOI = 25; // Trigger random normal generation every 25 iterations
                label currOIIt = 0; // Iteration block counter before random norm triggered

                for (label i=0; i<overlapIterations_; i++)
                {
                    // Delete the created molecule that is causing the overlap
                    if(newMol != NULL)
                    {
                        molCloud_.removeMolFromCellOccupancy(newMol->origId(), newMol->cell());
                        molCloud_.deleteParticle(*newMol);
                        newMol = NULL;
                    }

                    // Find perturbation for the position of the molecule based on the one it is overlapping
                    if(overlapMol != NULL)
                    {
                        overLap = position - overlapMol->position();
                        overLapMag = mag(overLap);

                        if(randomNorm == -1) //- Calculate normal if not using randomised values due to stuck molecule
                        {
                            overlapNorm = overLap / overLapMag;

                            //- Correct for rounding errors in normal calculation
                            if(overlapNorm[0] < -1.0)
                            {
                                overlapNorm[0] = -1.0;
                            }

                            if(overlapNorm[0] > 1.0)
                            {
                                overlapNorm[0] = 1.0;
                            }

                            if(overlapNorm[1] < -1.0)
                            {
                                overlapNorm[1] = -1.0;
                            }

                            if(overlapNorm[1] > 1.0)
                            {
                                overlapNorm[1] = 1.0;
                            }

                            if(overlapNorm[2] < -1.0)
                            {
                                overlapNorm[2] = -1.0;
                            }

                            if(overlapNorm[2] > 1.0)
                            {
                                overlapNorm[2] = 1.0;
                            }
                        }

                        if(overLapMag < 1.0)
                        {
                            perturbDistance = (1.0 / overLapMag) * perturbPerc; //Perturb molecule away from new overlapping molecule
                        }
                        else
                        {
                            perturbDistance = overLapMag * perturbPerc; //Perturb molecule away from new overlapping molecule
                        }

                        overlapMol = NULL;
                    }

                    vector origPos(position); //- Store the position before any perturbation

                    //Add perturbation to molecule position
                    position[0] += perturbDistance * overlapNorm[0];
                    position[1] += perturbDistance * overlapNorm[1];
                    position[2] += perturbDistance * overlapNorm[2];

                    //- If there are fixed boundary details then make sure the perturbed position doesn't fall outside of it
                    if(fixedBounds_)
                    {
                        trunkX = false;
                        trunkY = false;
                        trunkZ = false;

                        if(fixedBoundZeroThick_[0] == 1) //- Boundary has zero thickness in the x
                        {
                            if(fixedBoundNorm_[0] > 0)
                            {
                                if(position[0] <= fixedBoundMin_[0])
                                {
                                    position[0] = fixedBoundMin_[0] + boundCorr_;
                                    trunkX = true;
                                }
                            }
                            else if (fixedBoundNorm_[0] < 0)
                            {
                                if(position[0] >= fixedBoundMax_[0])
                                {
                                    position[0] = fixedBoundMax_[0] - boundCorr_;
                                    trunkX = true;
                                }
                            }
                        }
                        else
                        {
                            if(position[0] <= fixedBoundMin_[0])
                            {
                                position[0] = fixedBoundMin_[0] + boundCorr_;
                                trunkX = true;
                            }

                            if(position[0] >= fixedBoundMax_[0])
                            {
                                position[0] = fixedBoundMax_[0] - boundCorr_;
                                trunkX = true;
                            }
                        }

                        if(fixedBoundZeroThick_[1] == 1) //- Boundary has zero thickness in the y
                        {
                            if(fixedBoundNorm_[1] > 0)
                            {
                                if(position[1] <= fixedBoundMin_[1])
                                {
                                    position[1] = fixedBoundMin_[1] + boundCorr_;
                                    trunkY = true;
                                }
                            }
                            else if(fixedBoundNorm_[1] < 0)
                            {
                                if(position[1] >= fixedBoundMax_[1])
                                {
                                    position[1] = fixedBoundMax_[1] - boundCorr_;
                                    trunkY = true;
                                }
                            }
                        }
                        else
                        {
                            if(position[1] <= fixedBoundMin_[1])
                            {
                                position[1] = fixedBoundMin_[1] + boundCorr_;
                                trunkY = true;
                            }

                            if(position[1] >= fixedBoundMax_[1])
                            {
                                position[1] = fixedBoundMax_[1] - boundCorr_;
                                trunkY = true;
                            }
                        }

                        if(fixedBoundZeroThick_[2] == 1) //- Boundary has zero thickness in the z
                        {
                            if(fixedBoundNorm_[2] > 0)
                            {
                                if(position[2] <= fixedBoundMin_[2])
                                {
                                    position[2] = fixedBoundMin_[2] + boundCorr_;
                                    trunkZ = true;
                                }
                            }
                            else if(fixedBoundNorm_[2] < 0)
                            {
                                if(position[2] >= fixedBoundMax_[2])
                                {
                                    position[2] = fixedBoundMax_[2] - boundCorr_;
                                    trunkZ = true;
                                }
                            }
                        }
                        else
                        {
                            if(position[2] <= fixedBoundMin_[2])
                            {
                                position[2] = fixedBoundMin_[2] + boundCorr_;
                                trunkZ = true;
                            }

                            if(position[2] >= fixedBoundMax_[2])
                            {
                                position[2] = fixedBoundMax_[2] - boundCorr_;
                                trunkZ = true;
                            }
                        }

                        //- Check if truncation has happened in 2 or more directions (if so the particle is stuck)
                        if((trunkX && trunkY && trunkZ) ||
                           (trunkX && trunkY && !trunkZ) ||
                           (trunkX && !trunkY && trunkZ) ||
                           (!trunkX && trunkY && trunkZ))
                        {
                            stuckIter++;
                        }

                        //- Check if particle has been stuck for 2 iterations, if so reset its position and randomise normal to try again
                        if(stuckIter == 2)
                        {
                            randomNorm = 0; // Set randomNorm to 0 so a new random normal will be generated
                            stuckIter = 0; // Reset the stuckIter counter
                        }
                    }

                    if(randomNorm == 0) //- Calculate normal using random values
                    {
                        std::cout << "Generating random normal" << std::endl;

                        //- Undo the initial perturbation so new randomised normal can be used instead
                        position = origPos;

                        overlapNorm[0] = molCloud_.rndGen().position<scalar>(-1.0, 1.0);
                        overlapNorm[1] = molCloud_.rndGen().position<scalar>(-1.0, 1.0);
                        overlapNorm[2] = molCloud_.rndGen().position<scalar>(-1.0, 1.0);

                        //Add perturbation to molecule position using new randomised normal
                        position[0] += perturbDistance * overlapNorm[0];
                        position[1] += perturbDistance * overlapNorm[1];
                        position[2] += perturbDistance * overlapNorm[2];

                        //- If there are fixed boundary details then make sure the perturbed position doesn't fall outside of it
                        if(fixedBounds_)
                        {
                            if(fixedBoundZeroThick_[0] == 1) //- Boundary has zero thickness in the x
                            {
                                if(fixedBoundNorm_[0] > 0)
                                {
                                    if(position[0] <= fixedBoundMin_[0])
                                    {
                                        position[0] = fixedBoundMin_[0] + boundCorr_;
                                    }
                                }
                                else if (fixedBoundNorm_[0] < 0)
                                {
                                    if(position[0] >= fixedBoundMax_[0])
                                    {
                                        position[0] = fixedBoundMax_[0] - boundCorr_;
                                    }
                                }
                            }
                            else
                            {
                                if(position[0] <= fixedBoundMin_[0])
                                {
                                    position[0] = fixedBoundMin_[0] + boundCorr_;
                                }

                                if(position[0] >= fixedBoundMax_[0])
                                {
                                    position[0] = fixedBoundMax_[0] - boundCorr_;
                                }
                            }

                            if(fixedBoundZeroThick_[1] == 1) //- Boundary has zero thickness in the y
                            {
                                if(fixedBoundNorm_[1] > 0)
                                {
                                    if(position[1] <= fixedBoundMin_[1])
                                    {
                                        position[1] = fixedBoundMin_[1] + boundCorr_;
                                    }
                                }
                                else if(fixedBoundNorm_[1] < 0)
                                {
                                    if(position[1] >= fixedBoundMax_[1])
                                    {
                                        position[1] = fixedBoundMax_[1] - boundCorr_;
                                    }
                                }
                            }
                            else
                            {
                                if(position[1] <= fixedBoundMin_[1])
                                {
                                    position[1] = fixedBoundMin_[1] + boundCorr_;
                                }

                                if(position[1] >= fixedBoundMax_[1])
                                {
                                    position[1] = fixedBoundMax_[1] - boundCorr_;
                                }
                            }

                            if(fixedBoundZeroThick_[2] == 1) //- Boundary has zero thickness in the z
                            {
                                if(fixedBoundNorm_[2] > 0)
                                {
                                    if(position[2] <= fixedBoundMin_[2])
                                    {
                                        position[2] = fixedBoundMin_[2] + boundCorr_;
                                    }
                                }
                                else if(fixedBoundNorm_[2] < 0)
                                {
                                    if(position[2] >= fixedBoundMax_[2])
                                    {
                                        position[2] = fixedBoundMax_[2] - boundCorr_;
                                    }
                                }
                            }
                            else
                            {
                                if(position[2] <= fixedBoundMin_[2])
                                {
                                    position[2] = fixedBoundMin_[2] + boundCorr_;
                                }

                                if(position[2] >= fixedBoundMax_[2])
                                {
                                    position[2] = fixedBoundMax_[2] - boundCorr_;
                                }
                            }
                        }

                        //- Set randomNorm to 1 so a new random (or standard) norm won't be calculated (until a new stuck molecule is detected or iteration block reached)
                        randomNorm = 1;
                    }

                    //- Perturbation may have changed the cell, make sure the new position is a valid cell
                    cell = mesh_.findCell(position);

                    if(cell != -1)
                    {
                        tetFace = -1;
                        tetPt = -1;

                        // Update cell data before new insertion as perturbation may have moved things
                        mesh_.findCellFacePt
                        (
                            position,
                            cell,
                            tetFace,
                            tetPt
                        );

                        // Create in the newly perturbed location
                        newMol = molCloud_.createMolecule
                        (
                            position,
                            cell,
                            tetFace,
                            tetPt,
                            Q,
                            velocity,
                            vector::zero,
                            pi,
                            vector::zero,
                            specialPosition,
                            special,
                            id,
                            1.0,
                            molCloud_.getTrackingNumber()
                        );

                        molCloud_.insertMolInCellOccupancy(newMol);

                        overlapMol = checkForOverlaps(newMol, overlapEnergyLimit_);

                        iterCount++;

                        if(overlapMol == NULL) // The molecule no longer overlaps so break loop and carry on
                        {
                            break;
                        }
                    }
                    else
                    {
                        if(iterCount < overlapIterations_) //- Attempted cell was out of scope but there are more iterations remaining, reset position and trigger random normal
                        {
                            std::cout << "mdDsmcCoupling::insertMolecule(): Molecule insertion attempted outside of mesh whilst finding new location, trying again" << std::endl;
                            std::cout << "Attempted pos: " << position[0] << "," << position[1] << "," << position[2] << std::endl;
                            position = startPos;
                            randomNorm = 0;
                            iterCount++;
                        }
                        else //- Run out of iterations to find a new cell so abort insertion with warning
                        {
                            std::cout << "mdDsmcCoupling::insertMolecule(): Molecule insertion attempted outside of mesh whilst finding new location, molecule not inserted" << std::endl;
                            std::cout << "Attempted pos: " << position[0] << "," << position[1] << "," << position[2] << std::endl;
                            return NULL;
                        }
                    }

                    currOIIt++; //- Increment counter to trigger random norm calculation

                    if(currOIIt == randOI) //- Have reached the next block of the overall iteration count, trigger random norm generation
                    {
                        randomNorm = 0; //Set to zero so a new random norm will be generated during the next iteration
                        currOIIt = 0;
                    }
                }
            }

            if(overlapMol != NULL) //The iterative process to perturb the molecule away from overlap failed after nIter tries
            {
                if(newMol != NULL)
                {
                    // Delete the created molecule as it will exceed energy limit according to force-field calculation
                    molCloud_.removeMolFromCellOccupancy(newMol->origId(), newMol->cell());
                    molCloud_.deleteParticle(*newMol);
                    newMol = NULL;
                }

                std::cout << "mdDsmcCoupling::insertMolecule(): Failed to find new location for molecule. molecule not inserted" << std::endl;

                return NULL;
            }
            else if (iterCount > 0)
            {
                std::cout << "mdDsmcCoupling::insertMolecule(): New molecule insertion point found in " << iterCount << " iterations" << std::endl;
            }
        }

        return newMol;
    }
    else
    {
        std::cout << "mdDsmcCoupling::insertMolecule(): Molecule insertion attempted outside of mesh, molecule not inserted" << std::endl;
        std::cout << "Attempted pos: " << position[0] << "," << position[1] << "," << position[2] << std::endl;
        return NULL;
    }
}

polyMolecule* mdDsmcCoupling::checkForOverlaps(polyMolecule* newMol, const scalar& potEnergyLimit)
{
	polyMolecule* molJ = NULL;

	// Real-Real interactions

    forAll(molCloud_.il().dil(), d)
    {
        forAll(molCloud_.cellOccupancy()[d],cellIMols)
        {
            forAll(molCloud_.il().dil()[d], interactingCells)
            {
                List<polyMolecule*> cellJ =	molCloud_.cellOccupancy()[molCloud_.il().dil()[d][interactingCells]];

                forAll(cellJ, cellJMols)
                {
                    molJ = cellJ[cellJMols];

                     if(newMol->origId() != molJ->origId())
                     {
                         if(molCloud_.evaluatePotentialLimit(newMol, molJ, potEnergyLimit))
                         {
                             return molJ;
                         }
                     }
                }
            }
        }

        forAll(molCloud_.cellOccupancy()[d], cellIOtherMols)
        {
            molJ = molCloud_.cellOccupancy()[d][cellIOtherMols];

            if (molJ > newMol)
            {
                if(newMol->origId() != molJ->origId())
                {
                    if(molCloud_.evaluatePotentialLimit(newMol, molJ, potEnergyLimit))
                    {
                        return molJ;
                    }
                }
            }
        }
    }

    // Real-Referred interactions
    forAll(molCloud_.il().refCellsParticles(), r)
    {
        const List<label>& realCells = molCloud_.il().refCells()[r].neighbouringCells();

        forAll(molCloud_.il().refCellsParticles()[r], i)
        {
            molJ = molCloud_.il().refCellsParticles()[r][i];

            forAll(realCells, rC)
            {
                if(newMol->origId() != molJ->origId())
                {
                    if(molCloud_.evaluatePotentialLimit(newMol, molJ, potEnergyLimit))
                    {
                        return molJ;
                    }
                }
            }
        }
    }

	return NULL;
}

} // End namespace Foam

// ************************************************************************* //
