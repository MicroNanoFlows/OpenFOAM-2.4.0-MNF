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

#include "dsmcMdCoupling.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcMdCoupling, 0);

addToRunTimeSelectionTable(dsmcCouplingController, dsmcMdCoupling, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcMdCoupling::dsmcMdCoupling
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict,
    couplingInterface2d& twoDInterfaces,
    couplingInterface3d& threeDInterfaces
)
:
    dsmcCouplingController(t, cloud, dict, twoDInterfaces, threeDInterfaces),
    propsDict_(dict.subDict(typeName + "Properties")),
    propsDictSend_(dict.subDict(typeName + "Sending")),
    propsDictRecv_(dict.subDict(typeName + "Receiving")),
#ifdef USE_MUI
    cellCentres_(),
    sendInterfaces_(),
    recvInterfaces_(),
#endif
    sending_(false),
    receiving_(false),
    fixedBounds_(false),
    fixedRegion_(false),
    fixedRegionMin_(vector::zero),
    fixedRegionMax_(vector::zero),
    fixedBoundMin_(vector::zero),
    fixedBoundMax_(vector::zero),
    fixedBoundNorm_(vector::zero),
    fixedBoundZeroThick_(vector(-1, -1, -1)),
    boundCorr_(vector::zero),
    currIteration_(0)
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
        	Info << "dsmcMdCoupling::dsmcMdCoupling(): Found 3D MUI coupling interface ("
                 << interfaces[i] << ") to send for domain " << threeDInterfaces.domainName << endl;
        }

        parcelsInCellHistory_.setSize(sendInterfaces_.size());
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
        	Info << "dsmcMdCoupling::dsmcMdCoupling(): Found 3D MUI coupling interface ("
                 << interfaces[i] << ") to receive for domain " << threeDInterfaces.domainName << endl;
        }
    }

    if((sendInterfaces_.size() != 0))
    {
        sending_ = true;
    }

    if((recvInterfaces_.size() != 0))
    {
        receiving_ = true;
    }
#else
    FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                << "MUI library not enabled at compilation" << exit(FatalError);
#endif

    writeInTimeDir_ = true;
    writeInCase_ = true;

    const List<word> types(propsDict_.lookup("typeIds"));
    typeNames_.setSize(types.size());

    forAll(types, type)
    {
    	typeNames_[type] = types[type];
    }

    forAll(typeNames_, molType)
	{
    	const label typeId = findIndex(cloud_.typeIdList(), typeNames_[molType]);

		if(typeId == -1)
		{
			FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
						<< "Cannot find type id: " << typeNames_[molType] << nl << "in idList."
						<< exit(FatalError);
		}
	}

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
            FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                         << "A fixed boundary must have a zero thickness in at least one direction"
                         << exit(FatalError);
        }

        boundNormFound = true;
    }

    if((boundMinFound && !boundMaxFound && !boundNormFound) ||
       (boundMaxFound && !boundMinFound && !boundNormFound) ||
       (boundNormFound && !boundMinFound && !boundMaxFound))
    {
        FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                     << "Cannot find fixedBoundMin, fixedBoundMax and fixedBoundNorm together"
                     << exit(FatalError);
    }
    else
    {
        fixedBounds_ = true;

        const cellList& cells = mesh_.cells();
        const List<point>& pts = mesh_.points();

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

            if(fixedBoundZeroThick_[0] == 1)
            {
                if(fixedBoundNorm_[0] != 0)
                {
                    if(fixedBoundMin_[0] >= cellMin[0] && fixedBoundMax_[0] <= cellMax[0])
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
                    FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                                 << "Fixed boundary zero thickness in x direction but normal value zero"
                                 << exit(FatalError);
                }
            }

            if(fixedBoundZeroThick_[1] == 1)
            {
               if(fixedBoundNorm_[1] != 0)
               {
                   if(fixedBoundMin_[1] >= cellMin[1] && fixedBoundMax_[1] <= cellMax[1])
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
                   FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                                << "Fixed boundary zero thickness in y direction but normal value zero"
                                << exit(FatalError);
               }
            }

            if(fixedBoundZeroThick_[2] == 1)
            {
               if(fixedBoundNorm_[2] != 0)
               {
                   if(fixedBoundMin_[2] >= cellMin[2] && fixedBoundMax_[2] <= cellMax[2])
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
                   FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
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
                FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                             << "Fixed boundary defined but no intersecting cells found"
                             << exit(FatalError);
            }
        }
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcMdCoupling::~dsmcMdCoupling()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcMdCoupling::initialConfiguration()
{
    if(sending_)
    {
#ifdef USE_MUI
        boundCorr_[0] = ROOTVSMALL * 1e9;
        boundCorr_[1] = ROOTVSMALL * 1e9;
        boundCorr_[2] = ROOTVSMALL * 1e9;
#endif
        sendCoupledRegion(true); // Send ghost parcels in coupled regions at time = startTime
    }
}

void dsmcMdCoupling::controlParcelsBeforeCollisions()
{
    forAll(intersectingCells_, cell)
    {
        const List<dsmcParcel*>& parcelsInCell = cloud_.cellOccupancy()[intersectingCells_[cell]];

        forAll(parcelsInCell, parcel)
        {
            const word& parcType = cloud_.typeIdList()[parcelsInCell[parcel]->typeId()];
            label typeIndex = findIndex(typeNames_, parcType);

            //- Only delete and store particles that are of the coupled type
            if(typeIndex != -1)
            {
                bool removeParcel = false;

                if(fixedBoundZeroThick_[0] == 1)
                {
                    if(fixedBoundNorm_[0] > 0) //- Boundary is positive facing in the x
                    {
                        if(parcelsInCell[parcel]->position()[0] <= fixedBoundMin_[0])
                        {
                           removeParcel = true;
                        }
                    }
                    else if (fixedBoundNorm_[0] < 0) //- Boundary is negative facing in the x
                    {
                        if(parcelsInCell[parcel]->position()[0] >= fixedBoundMax_[0])
                        {
                            removeParcel = true;
                        }
                    }
                }

                if(fixedBoundZeroThick_[1] == 1)
                {
                    if(fixedBoundNorm_[1] > 0) //- Boundary is positive facing in the y
                    {
                        if(parcelsInCell[parcel]->position()[1] <= fixedBoundMin_[1])
                        {
                           removeParcel = true;
                        }
                    }
                    else if (fixedBoundNorm_[1] < 0) //- Boundary is negative facing in the x
                    {
                        if(parcelsInCell[parcel]->position()[1] >= fixedBoundMax_[1])
                        {
                            removeParcel = true;
                        }
                    }
                }

                if(fixedBoundZeroThick_[2] == 1)
                {
                    if(fixedBoundNorm_[2] > 0) //- Boundary is positive facing in the z
                    {
                        if(parcelsInCell[parcel]->position()[2] <= fixedBoundMin_[2])
                        {
                           removeParcel = true;
                        }
                    }
                    else if (fixedBoundNorm_[2] < 0) //- Boundary is negative facing in the x
                    {
                        if(parcelsInCell[parcel]->position()[2] >= fixedBoundMax_[2])
                        {
                            removeParcel = true;
                        }
                    }
                }

                if(removeParcel)
                {
                    //- Store the required details of the parcel
                    coupledParcel newParcToSend;
                    newParcToSend.parcType = cloud_.typeIdList()[parcelsInCell[parcel]->typeId()];
                    newParcToSend.position = parcelsInCell[parcel]->position();
                    newParcToSend.velocity = parcelsInCell[parcel]->U();
                    parcsToSend_.append(newParcToSend);

                    //- Delete parcel from cellOccupancy (before deleting it from cloud)
                    cloud_.removeParcelFromCellOccupancy(parcelsInCell[parcel]->origId(), intersectingCells_[cell]);
                    cloud_.deleteParticle(*parcelsInCell[parcel]);
                }
            }
        }
    }
}

void dsmcMdCoupling::calculateProperties(int stage)
{
    if(stage == 1)
    {
        currIteration_++; //- Increment current iteration

		if(sending_)
		{
			sendCoupledRegion(false); // Send ghost molecules in coupled regions (non-blocking)
		}
    }
    else if (stage == 2)
    {
    	if(receiving_)
		{
			receiveCoupledMolecules(); // Receive any molecules from MD coupling boundary (blocking)
		}
    }
    else if (stage == 3)
    {
    	if(sending_)
		{
			sendCoupledParcels(); // Send any parcels deleted by coupling boundary (non-blocking)
		}
    }
}

void dsmcMdCoupling::sendCoupledRegion(bool init)
{
#ifdef USE_MUI
    dsmcParcel* parcel = NULL;
    List<bool> listSizeChanged(sendInterfaces_.size(), false);
    bool parcelChanged = false;
    std::stringstream sendStr;

	// Iterate through all sending interfaces for this controller
	forAll(sendInterfaces_, iface)
	{
		DynamicList<label> currParcelsInCells;

		//First, create a current list of all parcels in cells for this interface to compare against historical list
		forAll(regionIds(), id)
		{
			forAll(controlZone(regionIds()[id]), c)
			{
				const label& cell = controlZone(regionIds()[id])[c];
				const List<dsmcParcel*>& parcelsInCell = cloud_.cellOccupancy()[cell];

				forAll(parcelsInCell, p) //Iterate through parcels in cell to determine if list size has changed
				{
					parcel = parcelsInCell[p];
                    currParcelsInCells.append(parcel->origId());
				}
			}
		}

		if(parcelsInCellHistory_[iface].size() != currParcelsInCells.size())
		{
			parcelsInCellHistory_[iface].clear();
			parcelsInCellHistory_[iface].setSize(currParcelsInCells.size());
			listSizeChanged[iface] = true;
		}

		forAll(regionIds(), id)
		{
			forAll(controlZone(regionIds()[id]), c)
			{
				const label& cell = controlZone(regionIds()[id])[c];
				const List<dsmcParcel*>& parcelsInCell = cloud_.cellOccupancy()[cell];

				forAll(parcelsInCell, p) // Iterate through parcels in cell
				{
					parcel = parcelsInCell[p];

					//- Determine whether parcel is of a type set to send
					const label typeIndex = findIndex(typeNames_, cloud_.typeIdList()[parcel->typeId()]);

					if(typeIndex != -1)
					{
						if(listSizeChanged[iface]) //List was resized so treat all parcels as new
						{
							parcelChanged = true;
							parcelsInCellHistory_[iface][p] = currParcelsInCells[p]; //Update history
						}
						else //List the same size as last time so only flag parcels that are different from the last time
						{
							if(parcelsInCellHistory_[iface][p] != currParcelsInCells[p])
							{
								parcelChanged = true;
								parcelsInCellHistory_[iface][p] = currParcelsInCells[p]; //Update history
							}
							else
							{
								parcelChanged = false;
							}
						}

						// Get the parcel centre
						mui::point3d parcCentre;
						parcCentre[0] = parcel->position()[0] * oneOverRefLength_;
						parcCentre[1] = parcel->position()[1] * oneOverRefLength_;
						parcCentre[2] = parcel->position()[2] * oneOverRefLength_;

                        // Push parcel type
                        sendInterfaces_[iface]->push("type_region", parcCentre, static_cast<std::string>(typeNames_[typeIndex]));

						// Push flag to say whether this parcel is new or not (to avoid creating molecules where they don't have to be)
						sendInterfaces_[iface]->push("changed_region", parcCentre, static_cast<scalar>(parcelChanged));

						// Push the parcel velocity to the interface
						sendInterfaces_[iface]->push("vel_x_region", parcCentre, parcel->U()[0]);
						sendInterfaces_[iface]->push("vel_y_region", parcCentre, parcel->U()[1]);
						sendInterfaces_[iface]->push("vel_z_region", parcCentre, parcel->U()[2]);
					}
				}
			}
		}

		std::cout << "Number of molecules sent to coupled region = " << currParcelsInCells.size() << std::endl;
	}

	forAll(sendInterfaces_, iface)
	{
		// Commit (transmit) values to the MUI interface
		sendInterfaces_[iface]->commit(currIteration_);
    }
#endif
}

void dsmcMdCoupling::sendCoupledParcels()
{
#ifdef USE_MUI
	std::stringstream sendStr;

	if(parcsToSend_.size() != 0)
	{
		forAll(parcsToSend_, parcs)
		{
            forAll(sendInterfaces_, iface)
            {
                // Get the parcel centre
                mui::point3d parcCentre;
                parcCentre[0] = parcsToSend_[parcs].position[0] * oneOverRefLength_;
                parcCentre[1] = parcsToSend_[parcs].position[1] * oneOverRefLength_;
                parcCentre[2] = parcsToSend_[parcs].position[2] * oneOverRefLength_;

                // Push parcel type
                sendInterfaces_[iface]->push("type_bound", parcCentre, static_cast<std::string>(parcsToSend_[parcs].parcType));

                // Push parcel velocity
                sendInterfaces_[iface]->push("vel_x_bound", parcCentre, parcsToSend_[parcs].velocity[0]);
                sendInterfaces_[iface]->push("vel_y_bound", parcCentre, parcsToSend_[parcs].velocity[1]);
                sendInterfaces_[iface]->push("vel_z_bound", parcCentre, parcsToSend_[parcs].velocity[2]);

                std::cout << "Coupling boundary parcel pushed at: " << "[" << parcCentre[0] << "," << parcCentre[1] << "," << parcCentre[2] << "]" << std::endl;
            }
		}

		//- Clear the sent parcels
		parcsToSend_.clear();
	}

	forAll(sendInterfaces_, iface)
	{
		// Commit (transmit) values to the MUI interface
		sendInterfaces_[iface]->commit(currIteration_);
	}
#endif
}

bool dsmcMdCoupling::receiveCoupledMolecules()
{
    bool parcelAdded = false;
#ifdef USE_MUI
	List<std::vector<mui::point3d> > rcvPoints(recvInterfaces_.size());
	List<std::vector<std::string> > rcvParcType(recvInterfaces_.size());
	List<std::vector<scalar> > rcvVelX(recvInterfaces_.size());
    List<std::vector<scalar> > rcvVelY(recvInterfaces_.size());
    List<std::vector<scalar> > rcvVelZ(recvInterfaces_.size());
	std::stringstream rcvStr;

	// Iterate through all receiving interfaces for this controller and extract a points list for each molecule type handled
    forAll(recvInterfaces_, iface)
    {
        //- Extract a list of all parcel locations
        rcvPoints[iface] = recvInterfaces_[iface]->fetch_points<std::string>("type_bound", currIteration_, *chrono_sampler);

        if(rcvPoints[iface].size() > 0)
        {
            //- Extract a list of all parcel types
            rcvParcType[iface] = recvInterfaces_[iface]->fetch_values<std::string>("type_bound", currIteration_, *chrono_sampler);

            //- Extract a list of all molecule velocities received from other solver through this interface
            rcvVelX[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_x_bound", currIteration_, *chrono_sampler);
            rcvVelY[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_y_bound", currIteration_, *chrono_sampler);
            rcvVelZ[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_z_bound", currIteration_, *chrono_sampler);
        }
    }

    List<DynamicList<label> > valuesToRemove(rcvParcType.size());

    //- Go through received values and find any that are not of the type set to be received
    forAll(rcvParcType, ifacepts)
    {
        if(rcvParcType[ifacepts].size() > 0)
        {
            for (size_t pts = 0; pts < rcvParcType[ifacepts].size(); ++pts)
            {
                const label typeIndex = findIndex(typeNames_, rcvParcType[ifacepts][pts]);

                if(typeIndex == -1) //- parcId not found in local list as one to receive so store it as one to remove from lists
                {
                    valuesToRemove[ifacepts].append(pts);
                }
            }
        }
    }

    //- Now remove any stored in the list not of a type to be received
    forAll(valuesToRemove, ifacepts)
    {
        if(valuesToRemove[ifacepts].size() > 0)
        {
            forAll(valuesToRemove[ifacepts], value)
            {
                rcvPoints[ifacepts].erase(rcvPoints[ifacepts].begin() + valuesToRemove[ifacepts][value]);
                rcvParcType[ifacepts].erase(rcvParcType[ifacepts].begin() + valuesToRemove[ifacepts][value]);
                rcvVelX[ifacepts].erase(rcvVelX[ifacepts].begin() + valuesToRemove[ifacepts][value]);
                rcvVelY[ifacepts].erase(rcvVelY[ifacepts].begin() + valuesToRemove[ifacepts][value]);
                rcvVelZ[ifacepts].erase(rcvVelZ[ifacepts].begin() + valuesToRemove[ifacepts][value]);
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
                velocity[0] = rcvVelX[ifacepts][pts];
                velocity[1] = rcvVelY[ifacepts][pts];
                velocity[2] = rcvVelZ[ifacepts][pts];

                point checkedPosition(rcvPoints[ifacepts][pts][0] * refLength_, rcvPoints[ifacepts][pts][1] * refLength_, rcvPoints[ifacepts][pts][2] * refLength_);

                if(fixedBounds_)
                {
                    if(fixedBoundZeroThick_[0] == 1) //- Boundary has zero thickness in the x
                    {
                        if(fixedBoundNorm_[0] != 0)
                        {
                            checkedPosition[0] = fixedBoundMin_[0] + (fixedBoundNorm_[0] * boundCorr_[0]); //- Move the new particle away from the boundary a little
                        }
                    }
                    else
                    {
                        if(checkedPosition[0] < fixedBoundMin_[0])
                        {
                            checkedPosition[0] = fixedBoundMin_[0] + boundCorr_[0];
                        }

                        if(checkedPosition[0] > fixedBoundMax_[0])
                        {
                            checkedPosition[0] = fixedBoundMax_[0] - boundCorr_[0];
                        }
                    }

                    if(fixedBoundZeroThick_[1] == 1) //- Boundary has zero thickness in the y
                    {
                        if(fixedBoundNorm_[1] != 0)
                        {
                            checkedPosition[1] = fixedBoundMin_[1] + (fixedBoundNorm_[1] * boundCorr_[1]); //- Move the new particle away from the boundary a little
                        }
                    }
                    else
                    {
                        if(checkedPosition[1] < fixedBoundMin_[1])
                        {
                            checkedPosition[1] = fixedBoundMin_[1] + boundCorr_[1];
                        }

                        if(checkedPosition[1] > fixedBoundMax_[1])
                        {
                            checkedPosition[1] = fixedBoundMax_[1] - boundCorr_[1];
                        }
                    }

                    if(fixedBoundZeroThick_[2] == 1) //- Boundary has zero thickness in the z
                    {
                        if(fixedBoundNorm_[2] != 0)
                        {
                            checkedPosition[2] = fixedBoundMin_[2] + (fixedBoundNorm_[2] * boundCorr_[2]); //- Move the new particle away from the boundary a little
                        }
                    }
                    else
                    {
                        if(checkedPosition[2] < fixedBoundMin_[2])
                        {
                            checkedPosition[2] = fixedBoundMin_[2] + boundCorr_[2];
                        }

                        if(checkedPosition[2] > fixedBoundMax_[2])
                        {
                            checkedPosition[2] = fixedBoundMax_[2] - boundCorr_[2];
                        }
                    }
                }

                const label typeIndex = findIndex(typeNames_, rcvParcType[ifacepts][pts]);

                insertParcel(checkedPosition, velocity, typeIndex);
                parcelAdded = true;

                std::cout << "Coupling boundary molecule inserted at: [" << checkedPosition[0] << "," << checkedPosition[1] << "," << checkedPosition[2] << "]" << std::endl;
            }
        }
	}
#endif
	return parcelAdded;
}

void dsmcMdCoupling::output
(
     const fileName& fixedPathName,
     const fileName& timePath
)
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        label singleStepNParcs = cloud_.size();
        scalar singleStepTotalLinearKE = 0.0;

        IDLList<dsmcParcel>::iterator parc(cloud_.begin());

        for(; parc != cloud_.end(); ++parc)
        {
            scalar parcMass(cloud_.constProps(parc().typeId()).mass());

            singleStepTotalLinearKE += 0.5*parcMass*magSqr(parc().U());
        }

        if (Pstream::parRun())
        {
            reduce(singleStepTotalLinearKE, sumOp<scalar>());
            reduce(singleStepNParcs, sumOp<label>());
        }

        if(Pstream::master())
        {
            OFstream os1(timePath/"avg_lin_KE");
            os1 << "Time " << time_.time().value() << endl;
            os1 << singleStepTotalLinearKE/singleStepNParcs << endl;
        }
    }
}

void dsmcMdCoupling::updateProperties(const dictionary& newDict)
{
	//- the main controller properties should be updated first
	updateCouplingControllerProperties(newDict);
	propsDict_ = newDict.subDict(typeName + "Properties");
	propsDictSend_ = newDict.subDict(typeName + "Sending");
	propsDictRecv_ = newDict.subDict(typeName + "Receiving");
}

void dsmcMdCoupling::barrier()
{
	forAll(sendInterfaces_, iface)
	{
		sendInterfaces_[iface]->barrier(currIteration_);
	}
}

void dsmcMdCoupling::barrier(label iteration)
{
	forAll(sendInterfaces_, iface)
	{
		sendInterfaces_[iface]->barrier(iteration);
	}
}

void dsmcMdCoupling::barrier(label iteration, label interface)
{
	// Wait for the other side to catch up
	sendInterfaces_[interface]->barrier(iteration);
}

void dsmcMdCoupling::forget()
{
	forAll(recvInterfaces_, iface)
	{
		recvInterfaces_[iface]->forget(currIteration_, true);
	}
}

void dsmcMdCoupling::forget(label iteration)
{
	forAll(recvInterfaces_, iface)
	{
		recvInterfaces_[iface]->forget(iteration, true);
	}
}

void dsmcMdCoupling::forget(label iteration, label interface)
{
	recvInterfaces_[interface]->forget(iteration, true);
}

void dsmcMdCoupling::insertParcel
(
    point& position,
	const vector& U,
	const label& typeId
)
{
    //These need to be properly defined, only correct for sample Argon case
    const scalar rotationalTemperature = 0;
    const scalar vibrationalTemperature = 0;
    const scalar electronicTemperature = 0;

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

        const dsmcParcel::constantProperties& cP = cloud_.constProps(typeId);

        scalar ERot = cloud_.equipartitionRotationalEnergy
        (
            rotationalTemperature,
            cP.rotationalDegreesOfFreedom()
        );

        labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
        (
            vibrationalTemperature,
            cP.vibrationalDegreesOfFreedom(),
            typeId
        );

        label ELevel = cloud_.equipartitionElectronicLevel
        (
            electronicTemperature,
            cP.degeneracyList(),
            cP.electronicEnergyList(),
            typeId
        );

        scalar RWF = 1.0;

        if(cloud_.axisymmetric())
        {
            const point& cC = mesh_.cellCentres()[cell];
            scalar radius = cC.y();

            RWF = 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());
        }

        scalarField wallTemperature(4, 0.0);
        vectorField wallVectors(4, vector::zero);

        cloud_.addNewParcel
        (
            position,
            U,
            RWF,
            ERot,
            ELevel,
            cell,
            tetFace,
            tetPt,
            typeId,
            0,
            0,
            0,
            wallTemperature,
            wallVectors,
            vibLevel 
        );
    }
    else
    {
        std::cout << "dsmcMdCoupling::insertParcel(): Parcel insertion attempted outside of mesh, parcel not inserted" << std::endl;
    }
}

} // End namespace Foam

// ************************************************************************* //
