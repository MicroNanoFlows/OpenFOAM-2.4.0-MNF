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
    sendingRegion_(false),
    receivingRegion_(false),
    sendingBound_(false),
    receivingBound_(false),
    couplingBounds_(false),
    couplingRegion_(false),
    couplingRegionMin_(vector::zero),
    couplingRegionMax_(vector::zero),
    couplingBoundMin_(vector::zero),
    couplingBoundMax_(vector::zero),
    couplingBoundNorm_(vector::zero),
    couplingBoundZeroThick_(vector(-1, -1, -1)),
    currIteration_(0),
    boundCorr_(0),
    meshMin_(VGREAT, VGREAT, VGREAT),
    meshMax_(-VSMALL, -VSMALL, -VSMALL)
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
        	std::cout << "dsmcMdCoupling::dsmcMdCoupling(): Found 3D MUI coupling interface ("
                 << interfaces[i] << ") to send for domain " << threeDInterfaces.domainName << std::endl;
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
        	std::cout << "dsmcMdCoupling::dsmcMdCoupling(): Found 3D MUI coupling interface ("
                 << interfaces[i] << ") to receive for domain " << threeDInterfaces.domainName << std::endl;
        }
    }

    if((sendInterfaces_.size() != 0))
    {
        sendingRegion_ = true;
        sendingBound_ = true;
    }

    if((recvInterfaces_.size() != 0))
    {
        receivingRegion_ = true;
        receivingBound_ = true;
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

    const List<point>& meshPoints = mesh_.points();

    //- Determine local mesh extents
    forAll(meshPoints, pts)
    {
        if(meshPoints[pts][0] < meshMin_[0])
        {
            meshMin_[0] = meshPoints[pts][0];
        }

        if(meshPoints[pts][1] < meshMin_[1])
        {
            meshMin_[1] = meshPoints[pts][1];
        }

        if(meshPoints[pts][2] < meshMin_[2])
        {
            meshMin_[2] = meshPoints[pts][2];
        }

        if(meshPoints[pts][0] > meshMax_[0])
        {
            meshMax_[0] = meshPoints[pts][0];
        }

        if(meshPoints[pts][1] > meshMax_[1])
        {
            meshMax_[1] = meshPoints[pts][1];
        }

        if(meshPoints[pts][2] > meshMax_[2])
        {
            meshMax_[2] = meshPoints[pts][2];
        }
    }

    bool regionMinFound = false;
    bool regionMaxFound = false;

    if (propsDict_.found("couplingRegionMin"))
    {
        couplingRegionMin_ = propsDict_.lookup("couplingRegionMin");
        regionMinFound = true;
    }

    if (propsDict_.found("couplingRegionMax"))
    {
        couplingRegionMax_ = propsDict_.lookup("couplingRegionMax");
        regionMaxFound = true;
    }

    const cellList& cells = mesh_.cells();

    if((regionMinFound && !regionMaxFound) || (regionMaxFound && !regionMinFound))
    {
        FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
                      << "Cannot find both couplingRegionMin and couplingRegionMax"
                      << exit(FatalError);
    }
    else
    {
        couplingRegion_ = true;

        vector meshHalfWidth(((meshMax_[0] - meshMin_[0]) * 0.5),
                             ((meshMax_[1] - meshMin_[1]) * 0.5),
                             ((meshMax_[2] - meshMin_[2]) * 0.5));
        vector couplingRegionHalfWidth(((couplingRegionMax_[0] - couplingRegionMin_[0]) * 0.5),
                                   ((couplingRegionMax_[1] - couplingRegionMin_[1]) * 0.5),
                                   ((couplingRegionMax_[2] - couplingRegionMin_[2]) * 0.5));
        point meshCentre(meshMin_[0] + meshHalfWidth[0],
                         meshMin_[1] + meshHalfWidth[1],
                         meshMin_[2] + meshHalfWidth[2]);
        point couplingRegionCentre(couplingRegionMin_[0] + couplingRegionHalfWidth[0],
                                   couplingRegionMin_[1] + couplingRegionHalfWidth[1],
                                   couplingRegionMin_[2] + couplingRegionHalfWidth[2]);

        bool overlap = true;

        if ((std::fabs(meshCentre[0] - couplingRegionCentre[0]) > (meshHalfWidth[0] + couplingRegionHalfWidth[0])) ||
           (std::fabs(meshCentre[1] - couplingRegionCentre[1]) > (meshHalfWidth[1] + couplingRegionHalfWidth[1])) ||
           (std::fabs(meshCentre[2] - couplingRegionCentre[2]) > (meshHalfWidth[2] + couplingRegionHalfWidth[2])))
        {
            overlap = false;
        }

        //- There is an overlap between the coupling region and the local mesh so should have found at least 1 intersecting cell
        if(!overlap)
        {
            sendingRegion_ = false;
            receivingRegion_ = false;
        }
        else //- Found an overlap between the region and local mesh so find which cells the region is in
        {
            point cellMin;
            point cellMax;

            //- Determine which cells the coupling region intersects
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
                    if(meshPoints[pointList[cellPoint]][0] < cellMin[0])
                    {
                        cellMin[0] = meshPoints[pointList[cellPoint]][0];
                    }

                    if(meshPoints[pointList[cellPoint]][0] > cellMax[0])
                    {
                        cellMax[0] = meshPoints[pointList[cellPoint]][0];
                    }

                    if(meshPoints[pointList[cellPoint]][1] < cellMin[1])
                    {
                        cellMin[1] = meshPoints[pointList[cellPoint]][1];
                    }

                    if(meshPoints[pointList[cellPoint]][1] > cellMax[1])
                    {
                        cellMax[1] = meshPoints[pointList[cellPoint]][1];
                    }

                    if(meshPoints[pointList[cellPoint]][2] < cellMin[2])
                    {
                        cellMin[2] = meshPoints[pointList[cellPoint]][2];
                    }

                    if(meshPoints[pointList[cellPoint]][2] > cellMax[2])
                    {
                        cellMax[2] = meshPoints[pointList[cellPoint]][2];
                    }
                }

                vector cellHalfWidth(((cellMax[0] - cellMin[0]) * 0.5),
                                     ((cellMax[1] - cellMin[1]) * 0.5),
                                     ((cellMax[2] - cellMin[2]) * 0.5));
                point cellCentre(cellMin[0] + cellHalfWidth[0],
                                 cellMin[1] + cellHalfWidth[1],
                                 cellMin[2] + cellHalfWidth[2]);

                bool overlap = true;

                if ((std::fabs(cellCentre[0] - couplingRegionCentre[0]) > (cellHalfWidth[0] + couplingRegionHalfWidth[0])) ||
                   (std::fabs(cellCentre[1] - couplingRegionCentre[1]) > (cellHalfWidth[1] + couplingRegionHalfWidth[1])) ||
                   (std::fabs(cellCentre[2] - couplingRegionCentre[2]) > (cellHalfWidth[2] + couplingRegionHalfWidth[2])))
                {
                    overlap = false;
                }

                if(overlap)
                {
                    regionCells_.append(cell);
                }
            }

            if(regionCells_.size() > 0)
            {
                std::cout << "Region intersecting cell count: " << regionCells_.size() << std::endl;
            }
        }
    }

    bool boundMinFound = false;
    bool boundMaxFound = false;
    bool boundNormFound = false;

    if (propsDict_.found("couplingBoundMin"))
    {
        couplingBoundMin_ = propsDict_.lookup("couplingBoundMin");
        boundMinFound = true;
    }

    if (propsDict_.found("couplingBoundMax"))
    {
        couplingBoundMax_ = propsDict_.lookup("couplingBoundMax");
        boundMaxFound = true;
    }

    if (propsDict_.found("couplingBoundNorm"))
    {
        couplingBoundNorm_ = propsDict_.lookup("couplingBoundNorm");

        if(couplingBoundMax_[0] - couplingBoundMin_[0] == 0.0)
        {
            couplingBoundZeroThick_[0] = 1;
        }

        if(couplingBoundMax_[1] - couplingBoundMin_[1] == 0.0)
        {
            couplingBoundZeroThick_[1] = 1;
        }

        if(couplingBoundMax_[2] - couplingBoundMin_[2] == 0.0)
        {
            couplingBoundZeroThick_[2] = 1;
        }

        if(couplingBoundZeroThick_[0] == -1 && couplingBoundZeroThick_[1] == -1 && couplingBoundZeroThick_[2] == -1)
        {
            FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                         << "A coupling boundary must have a zero thickness in at least one direction"
                         << exit(FatalError);
        }

        boundNormFound = true;
    }

    if((boundMinFound && !boundMaxFound && !boundNormFound) ||
       (boundMaxFound && !boundMinFound && !boundNormFound) ||
       (boundNormFound && !boundMinFound && !boundMaxFound))
    {
        FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                     << "Cannot find couplingBoundMin, couplingBoundMax and couplingBoundNorm together"
                     << exit(FatalError);
    }
    else
    {
        couplingBounds_ = true;

        vector localMeshExtents = meshMax_ - meshMin_;

        vector meshExtents = mesh_.bounds().max() - mesh_.bounds().min();

        vector boundCorr = meshExtents * 1e-8;

        //- Ensure boundary correction value not larger than 1e-9
        if(boundCorr[0] > 1e-9)
        {
            boundCorr[0] = 1e-9;
        }

        if(boundCorr[1] > 1e-9)
        {
            boundCorr[1] = 1e-9;
        }

        if(boundCorr[2] > 1e-9)
        {
            boundCorr[2] = 1e-9;
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
        bool intersectCell = false;

        //- Determine which cells the coupling boundary intersects
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
                if(meshPoints[pointList[cellPoint]][0] < cellMin[0])
                {
                    cellMin[0] = meshPoints[pointList[cellPoint]][0];
                }

                if(meshPoints[pointList[cellPoint]][0] > cellMax[0])
                {
                    cellMax[0] = meshPoints[pointList[cellPoint]][0];
                }

                if(meshPoints[pointList[cellPoint]][1] < cellMin[1])
                {
                    cellMin[1] = meshPoints[pointList[cellPoint]][1];
                }

                if(meshPoints[pointList[cellPoint]][1] > cellMax[1])
                {
                    cellMax[1] = meshPoints[pointList[cellPoint]][1];
                }

                if(meshPoints[pointList[cellPoint]][2] < cellMin[2])
                {
                    cellMin[2] = meshPoints[pointList[cellPoint]][2];
                }

                if(meshPoints[pointList[cellPoint]][2] > cellMax[2])
                {
                    cellMax[2] = meshPoints[pointList[cellPoint]][2];
                }
            }

            if(couplingBoundZeroThick_[0] == 1) //- couplingBoundMin_ and couplingBoundMax_ are the same in the x
            {
                if(couplingBoundNorm_[0] != 0)
                {
                    scalar boundaryExtend = localMeshExtents[0] * 1e-4;

                    bool test = false;

                    if((couplingBoundMin_[0] - boundaryExtend) >= cellMin[0] && (couplingBoundMin_[0] - boundaryExtend) <= cellMax[0])
                    {
                        test = true;
                    }

                    if((couplingBoundMin_[0] + boundaryExtend) >= cellMin[0] && (couplingBoundMin_[0] + boundaryExtend) <= cellMax[0])
                    {
                        test = true;
                    }

                    if(test)
                    {
                        if((cellMin[1] >= couplingBoundMin_[1] && cellMax[1] <= couplingBoundMax_[1]) &&
                           (cellMin[2] >= couplingBoundMin_[2] && cellMax[2] <= couplingBoundMax_[2]))
                        {
                            intersectCell = true;
                        }
                    }
                }
                else
                {
                    FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                                 << "Coupling boundary zero thickness in x direction but normal value zero"
                                 << exit(FatalError);
                }
            }

            if(couplingBoundZeroThick_[1] == 1) //- couplingBoundMin_ and couplingBoundMax_ are the same in the y
            {
               if(couplingBoundNorm_[1] != 0)
               {
                   scalar boundaryExtend = localMeshExtents[1] * 1e-4;
                   bool test = false;

                   if((couplingBoundMin_[1] - boundaryExtend) >= cellMin[1] && (couplingBoundMin_[1] - boundaryExtend) <= cellMax[1])
                   {
                       test = true;
                   }

                   if((couplingBoundMin_[1] + boundaryExtend) >= cellMin[1] && (couplingBoundMin_[1] + boundaryExtend) <= cellMax[1])
                   {
                       test = true;
                   }

                   if(test)
                   {
                       if((cellMin[0] >= couplingBoundMin_[0] && cellMax[0] <= couplingBoundMax_[0]) &&
                          (cellMin[2] >= couplingBoundMin_[2] && cellMax[2] <= couplingBoundMax_[2]))
                       {
                           intersectCell = true;
                       }
                   }
               }
               else
               {
                   FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                                << "Coupling boundary zero thickness in y direction but normal value zero"
                                << exit(FatalError);
               }
            }

            if(couplingBoundZeroThick_[2] == 1) //- couplingBoundMin_ and couplingBoundMax_ are the same in the z
            {
               if(couplingBoundNorm_[2] != 0)
               {
                   scalar boundaryExtend = localMeshExtents[2] * 1e-4;
                   bool test = false;

                   if((couplingBoundMin_[2] - boundaryExtend) >= cellMin[2] && (couplingBoundMin_[2] - boundaryExtend) <= cellMax[2])
                   {
                       test = true;
                   }

                   if((couplingBoundMin_[2] + boundaryExtend) >= cellMin[2] && (couplingBoundMin_[2] + boundaryExtend) <= cellMax[2])
                   {
                       test = true;
                   }

                   if(test)
                   {
                       if((cellMin[0] >= couplingBoundMin_[0] && cellMax[0] <= couplingBoundMax_[0]) &&
                          (cellMin[1] >= couplingBoundMin_[1] && cellMax[1] <= couplingBoundMax_[1]))
                       {
                           intersectCell = true;
                       }
                   }
               }
               else
               {
                   FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                                << "Coupling boundary zero thickness in z direction but normal value zero"
                                << exit(FatalError);
               }
            }
        }

        if(!intersectCell)
        {
            vector meshHalfWidth(((meshMax_[0] - meshMin_[0]) * 0.5),
                                 ((meshMax_[1] - meshMin_[1]) * 0.5),
                                 ((meshMax_[2] - meshMin_[2]) * 0.5));
            vector couplingBoundHalfWidth(((couplingBoundMax_[0] - couplingBoundMin_[0]) * 0.5),
                                       ((couplingBoundMax_[1] - couplingBoundMin_[1]) * 0.5),
                                       ((couplingBoundMax_[2] - couplingBoundMin_[2]) * 0.5));
            point meshCentre(meshMin_[0] + meshHalfWidth[0],
                             meshMin_[1] + meshHalfWidth[1],
                             meshMin_[2] + meshHalfWidth[2]);
            point couplingBoundCentre(couplingBoundMin_[0] + couplingBoundHalfWidth[0],
                                   couplingBoundMin_[1] + couplingBoundHalfWidth[1],
                                   couplingBoundMin_[2] + couplingBoundHalfWidth[2]);

            bool overlap = true;

            if ((std::fabs(meshCentre[0] - couplingBoundCentre[0]) > (meshHalfWidth[0] + couplingBoundHalfWidth[0])) ||
               (std::fabs(meshCentre[1] - couplingBoundCentre[1]) > (meshHalfWidth[1] + couplingBoundHalfWidth[1])) ||
               (std::fabs(meshCentre[2] - couplingBoundCentre[2]) > (meshHalfWidth[2] + couplingBoundHalfWidth[2])))
            {
                overlap = false;
            }

            //- There is an overlap between the coupling boundary and the local mesh so should have found at least 1 intersecting cell
            if(overlap)
            {
                FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                             << "Coupling boundary defined but no intersecting cells found"
                             << exit(FatalError);
            }
            else
            {
                sendingBound_ = false;
                receivingBound_ = false;
            }
        }
    }

    refLength_ = threeDInterfaces.refLength; //- Store the reference length
    refTime_ = threeDInterfaces.refTime; //- Store the reference time

    oneOverRefLength_ = 1.0 / refLength_;
    oneOverRefTime_ = 1.0 / refTime_;
#ifdef USE_MUI
    //Initialise exact time sampler for MUI
    chrono_sampler = new mui::chrono_sampler_exact3d();
#endif
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcMdCoupling::~dsmcMdCoupling()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcMdCoupling::initialConfiguration(label stage)
{
    if(stage == 1)
    {
#ifdef USE_MUI
    if((!sendingBound_ && !sendingRegion_) && (!receivingBound_ && !receivingRegion_))
    {
        std::cout << "MUI interface(s) disabled for this rank" << std::endl;

        forAll(sendInterfaces_, iface)
        {
            sendInterfaces_[iface]->announce_send_disable();
        }

        forAll(recvInterfaces_, iface)
        {
            recvInterfaces_[iface]->announce_recv_disable();
        }
    }
    else if((!sendingBound_ && !sendingRegion_))
    {
        std::cout << "MUI interface(s) sending disabled for this rank" << std::endl;

        forAll(sendInterfaces_, iface)
        {
            sendInterfaces_[iface]->announce_send_disable();
        }
    }
    else if((!receivingBound_ && !receivingRegion_))
    {
        std::cout << "MUI interface(s) receiving disabled for this rank" << std::endl;

        forAll(recvInterfaces_, iface)
        {
            recvInterfaces_[iface]->announce_recv_disable();
        }
    }

    DynamicList<word> interfaceCommits;

    forAll(sendInterfaces_, iface)
    {
        sendInterfaces_[iface]->commit(static_cast<label>(-1));
        interfaceCommits.append(sendInterfaceNames_[iface]);
    }

    forAll(recvInterfaces_, iface)
    {
        label index = findIndex(interfaceCommits, recvInterfaceNames_[iface]);
        if(index == -1)
        {
            recvInterfaces_[iface]->commit(static_cast<label>(-1));
        }
    }
#endif
    }
    else if (stage == 2)
    {
        sendCoupledRegion(true); // Send ghost parcels in coupled regions at time = startTime
    }
}

void dsmcMdCoupling::controlParcelsBeforeCollisions(label stage)
{
    if(stage == 1) //- Collate and delete parcels that have passed through a coupling boundary
    {
        currIteration_++; //- Increment current iteration

        if(receiveCoupledMolecules()) // Receive any molecules from MD coupling boundary (blocking)
        {
            // Update cell occupancy if parcels were received and inserted
            cloud_.reBuildCellOccupancy();
        }
    }
    else if (stage == 2) //- Send parcels that were deleted by coupling boundary during evolve() in dsmcCloud
    {
        sendCoupledParcels(); // Send any parcels deleted by coupling boundary (non-blocking)
    }
    else if (stage == 3) //- Clear the sent parcels in the cloud
    {
        cloud_.clearCoupledParcels();
    }
}

void dsmcMdCoupling::controlParcelsAfterCollisions(int stage)
{
    if(stage == 1)
    {
        sendCoupledRegion(false); // Send parcel positions in coupled region(s) (non-blocking)
	}
    else if (stage == 2)
    {
        receiveCoupledRegionForces(); // Receive MD forces on ghost molecules in coupled region(s) (blocking)
	}
}

void dsmcMdCoupling::sendCoupledRegion(bool init)
{
#ifdef USE_MUI
    dsmcParcel* parcel = NULL;

	// Iterate through all sending interfaces for this controller
    if(sendingRegion_)
    {
        forAll(sendInterfaces_, iface)
        {
            parcelsInCellHistory_[iface].clear(); // Clear the send history list

            forAll(regionCells_, cell)
            {
                const List<dsmcParcel*>& parcelsInCell = cloud_.cellOccupancy()[regionCells_[cell]];

                forAll(parcelsInCell, p)
                {
                    parcel = parcelsInCell[p];
                    bool insideRegion = false;

                    if(couplingRegion_)
                    {
                        vector position = parcel->position();

                        if((position[0] > couplingRegionMin_[0] && position[0] < couplingRegionMax_[0]) &&
                           (position[1] > couplingRegionMin_[1] && position[1] < couplingRegionMax_[1]) &&
                           (position[2] > couplingRegionMin_[2] && position[2] < couplingRegionMax_[2]))
                        {
                            insideRegion = true;
                        }
                    }
                    else
                    {
                        insideRegion = true;
                    }

                    if(insideRegion)
                    {
                        //- Determine whether parcel is of a type set to send
                        const label typeIndex = findIndex(typeNames_, cloud_.typeIdList()[parcel->typeId()]);

                        if(typeIndex != -1)
                        {
                            // Get the parcel centre
                            mui::point3d parcCentre;
                            parcCentre[0] = parcel->position()[0] * oneOverRefLength_;
                            parcCentre[1] = parcel->position()[1] * oneOverRefLength_;
                            parcCentre[2] = parcel->position()[2] * oneOverRefLength_;

                            // Push parcel type
                            sendInterfaces_[iface]->push("type_region", parcCentre, static_cast<std::string>(typeNames_[typeIndex]));

                            // Push parcel ID
                            sendInterfaces_[iface]->push("id_region", parcCentre, static_cast<label>(parcel->origId()));

                            // Push the parcel velocity to the interface
                            sendInterfaces_[iface]->push("vel_x_region", parcCentre, parcel->U()[0]);
                            sendInterfaces_[iface]->push("vel_y_region", parcCentre, parcel->U()[1]);
                            sendInterfaces_[iface]->push("vel_z_region", parcCentre, parcel->U()[2]);

                            parcelsInCellHistory_[iface].append(parcel->origId()); //Store the parcel ID that was sent
                        }
                    }
                }
            }

            std::cout << "    Parcels sent to coupled region  = " << parcelsInCellHistory_[iface].size() << std::endl;

            // Commit (transmit) values to the MUI interface
            sendInterfaces_[iface]->commit(currIteration_);
	    }
	}
#endif
}

void dsmcMdCoupling::receiveCoupledRegionForces()
{
#ifdef USE_MUI
    if(receivingRegion_)
    {
        List<std::vector<std::string> > rcvParcType(recvInterfaces_.size());
        List<std::vector<label> > rcvParcId(recvInterfaces_.size());
        List<std::vector<scalar> > rcvForceX(recvInterfaces_.size());
        List<std::vector<scalar> > rcvForceY(recvInterfaces_.size());
        List<std::vector<scalar> > rcvForceZ(recvInterfaces_.size());

        // Iterate through all receiving interfaces for this controller and extract a points list
        forAll(recvInterfaces_, iface)
        {
            //- Extract a list of all parcel types
            rcvParcType[iface] = recvInterfaces_[iface]->fetch_values<std::string>("type_region", currIteration_, *chrono_sampler);

            if(rcvParcType[iface].size() > 0)
            {
                //- Extract a list of all molecule Id's received from other solver through this interface
                rcvParcId[iface] = recvInterfaces_[iface]->fetch_values<label>("id_region", currIteration_, *chrono_sampler);

                //- Extract a list of all molecule forces received from other solver through this interface
                rcvForceX[iface] = recvInterfaces_[iface]->fetch_values<scalar>("force_x_region", currIteration_, *chrono_sampler);
                rcvForceY[iface] = recvInterfaces_[iface]->fetch_values<scalar>("force_y_region", currIteration_, *chrono_sampler);
                rcvForceZ[iface] = recvInterfaces_[iface]->fetch_values<scalar>("force_z_region", currIteration_, *chrono_sampler);
            }
        }

        //- Go through received values and find any that are not of the type set to be received
        forAll(rcvParcType, ifacepts)
        {
            if(rcvParcType[ifacepts].size() > 0)
            {
                std::vector<std::string>::iterator rcvParcTypeIt;
                std::vector<label>::iterator rcvParcIdIt = rcvParcId[ifacepts].begin();
                std::vector<scalar>::iterator rcvForceXIt = rcvForceX[ifacepts].begin();
                std::vector<scalar>::iterator rcvForceYIt = rcvForceY[ifacepts].begin();
                std::vector<scalar>::iterator rcvForceZIt = rcvForceZ[ifacepts].begin();

                for (rcvParcTypeIt = rcvParcType[ifacepts].begin(); rcvParcTypeIt != rcvParcType[ifacepts].end(); rcvParcTypeIt++) {
                    const label parcId = findIndex(typeNames_, *rcvParcTypeIt);

                    if(parcId == -1) //- parcId not found in local list as one to receive so store it as one to remove from lists
                    {
                        rcvParcType[ifacepts].erase(rcvParcTypeIt--);
                        rcvParcId[ifacepts].erase(rcvParcIdIt--);
                        rcvForceX[ifacepts].erase(rcvForceXIt--);
                        rcvForceY[ifacepts].erase(rcvForceYIt--);
                        rcvForceZ[ifacepts].erase(rcvForceZIt--);
                    }

                    rcvParcIdIt++;
                    rcvForceXIt++;
                    rcvForceYIt++;
                    rcvForceZIt++;
                }
            }
        }

        DynamicList<dsmcParcel*> parcelsInRegion;
        dsmcParcel* parcel = NULL;

        //- Extract a list of all parcels in the control zone(s)
        forAll(regionCells_, cell)
        {
            const List<dsmcParcel*>& parcelsInCell = cloud_.cellOccupancy()[regionCells_[cell]];

            forAll(parcelsInCell, p)
            {
                parcel = parcelsInCell[p];
                bool insideRegion = false;

                if(couplingRegion_)
                {
                    vector position = parcel->position();

                    if((position[0] > couplingRegionMin_[0] && position[0] < couplingRegionMax_[0]) &&
                       (position[1] > couplingRegionMin_[1] && position[1] < couplingRegionMax_[1]) &&
                       (position[2] > couplingRegionMin_[2] && position[2] < couplingRegionMax_[2]))
                    {
                        insideRegion = true;
                    }
                }
                else
                {
                    insideRegion = true;
                }

                if(insideRegion)
                {
                    //- Determine whether parcel is of a type set to receive
                    const label typeIndex = findIndex(typeNames_, cloud_.typeIdList()[parcel->typeId()]);

                    if(typeIndex != -1)
                    {
                        parcelsInRegion.append(parcel);
                    }
                }
            }
        }


        // Iterate through all forces received for this controller and apply if IDs match
        forAll(rcvParcId, iface)
        {
            if(parcelsInRegion.size() == rcvParcId[iface].size())
            {
                forAll(parcelsInRegion, parcel)
                {
                    if(rcvParcId[iface][parcel] == parcelsInRegion[parcel]->origId())
                    {
                        scalar parcMass(cloud_.constProps(parcelsInRegion[parcel]->typeId()).mass());

                        parcelsInRegion[parcel]->U()[0] += 0.5 * (rcvForceX[iface][parcel] / parcMass) * mesh_.time().deltaTValue();
                        parcelsInRegion[parcel]->U()[1] += 0.5 * (rcvForceY[iface][parcel] / parcMass) * mesh_.time().deltaTValue();
                        parcelsInRegion[parcel]->U()[2] += 0.5 * (rcvForceZ[iface][parcel] / parcMass) * mesh_.time().deltaTValue();
                    }
                }
            }
        }
    }
#endif
}

void dsmcMdCoupling::sendCoupledParcels()
{
#ifdef USE_MUI
    if(sendingBound_)
    {
        dsmcParcel* parc = NULL;

        forAll(sendInterfaces_, iface)
        {
            label pushed = 0;
            const DynamicList<dsmcCloud::coupledParc>& parcsToSend = cloud_.coupledParcels();

            if(parcsToSend.size() > 0)
            {
                forAll(parcsToSend, parcs)
                {
                    //- Only send the parcel if the sending interface defined in the boundary matches current interface
                    if(findIndex(parcsToSend[parcs].sendingInterfaces, sendInterfaceNames_[iface]) != -1)
                    {
                        parc = parcsToSend[parcs].parcel;

                        // Get the parcel centre
                        mui::point3d parcCentre;
                        parcCentre[0] = parc->position()[0] * oneOverRefLength_;
                        parcCentre[1] = parc->position()[1] * oneOverRefLength_;
                        parcCentre[2] = parc->position()[2] * oneOverRefLength_;

                        // Push parcel type
                        sendInterfaces_[iface]->push("type_bound", parcCentre, static_cast<std::string>(cloud_.typeIdList()[parc->typeId()]));

                        // Push parcel velocity
                        sendInterfaces_[iface]->push("vel_x_bound", parcCentre, parc->U()[0]);
                        sendInterfaces_[iface]->push("vel_y_bound", parcCentre, parc->U()[1]);
                        sendInterfaces_[iface]->push("vel_z_bound", parcCentre, parc->U()[2]);

                        pushed++;
                    }
                }
            }

            if(pushed > 0)
            {
                std::cout << "    Coupling parcels pushed         = " << pushed << std::endl;
            }

            // Commit (transmit) values to the MUI interface
            sendInterfaces_[iface]->commit(currIteration_);
        }
	}
#endif
}

bool dsmcMdCoupling::receiveCoupledMolecules()
{
    bool parcelAdded = false;
#ifdef USE_MUI
    if(receivingBound_)
    {
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

        //- Go through received values and find any that are not of the type set to be received
        forAll(rcvParcType, ifacepts)
        {
            if(rcvParcType[ifacepts].size() > 0)
            {
                std::vector<std::string>::iterator rcvParcTypeIt;
                std::vector<mui::point3d>::iterator rcvPointsIt = rcvPoints[ifacepts].begin();
                std::vector<scalar>::iterator rcvVelXIt = rcvVelX[ifacepts].begin();
                std::vector<scalar>::iterator rcvVelYIt = rcvVelY[ifacepts].begin();
                std::vector<scalar>::iterator rcvVelZIt = rcvVelZ[ifacepts].begin();

                for (rcvParcTypeIt = rcvParcType[ifacepts].begin(); rcvParcTypeIt != rcvParcType[ifacepts].end(); rcvParcTypeIt++) {
                    const label parcId = findIndex(typeNames_, *rcvParcTypeIt);

                    if(parcId == -1) //- parcId not found in local list as one to receive so store it as one to remove from lists
                    {
                        rcvParcType[ifacepts].erase(rcvParcTypeIt--);
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

        label inserted = 0;

        // Iterate through all received points
        forAll(rcvPoints, ifacepts)
        {
            if(rcvPoints[ifacepts].size() > 0)
            {
                for (size_t pts = 0; pts < rcvPoints[ifacepts].size(); pts++)
                {
                    point checkedPosition(rcvPoints[ifacepts][pts][0] * refLength_, rcvPoints[ifacepts][pts][1] * refLength_, rcvPoints[ifacepts][pts][2] * refLength_);

                    label cell = mesh_.findCell(checkedPosition);

                    // Attempt insertion/update only if received molecule falls within local mesh extents
                    if(cell != -1)
                    {
                        vector velocity;
                        velocity[0] = rcvVelX[ifacepts][pts];
                        velocity[1] = rcvVelY[ifacepts][pts];
                        velocity[2] = rcvVelZ[ifacepts][pts];

                        if(couplingBounds_)
                        {
                            if(couplingBoundZeroThick_[0] == 1) //- Boundary has zero thickness in the x
                            {
                                checkedPosition[0] = couplingBoundMin_[0] + (couplingBoundNorm_[0] * boundCorr_);
                            }
                            else
                            {
                                if(checkedPosition[0] <= couplingBoundMin_[0])
                                {
                                    checkedPosition[0] = couplingBoundMin_[0] + boundCorr_;
                                }

                                if(checkedPosition[0] >= couplingBoundMax_[0])
                                {
                                    checkedPosition[0] = couplingBoundMax_[0] - boundCorr_;
                                }
                            }

                            if(couplingBoundZeroThick_[1] == 1) //- Boundary has zero thickness in the y
                            {
                                checkedPosition[1] = couplingBoundMin_[1] + (couplingBoundNorm_[1] * boundCorr_);
                            }
                            else
                            {
                                if(checkedPosition[1] <= couplingBoundMin_[1])
                                {
                                    checkedPosition[1] = couplingBoundMin_[1] + boundCorr_;
                                }

                                if(checkedPosition[1] >= couplingBoundMax_[1])
                                {
                                    checkedPosition[1] = couplingBoundMax_[1] - boundCorr_;
                                }
                            }

                            if(couplingBoundZeroThick_[2] == 1) //- Boundary has zero thickness in the z
                            {
                                checkedPosition[2] = couplingBoundMin_[2] + (couplingBoundNorm_[2] * boundCorr_);
                            }
                            else
                            {
                                if(checkedPosition[2] <= couplingBoundMin_[2])
                                {
                                    checkedPosition[2] = couplingBoundMin_[2] + boundCorr_;
                                }

                                if(checkedPosition[2] >= couplingBoundMax_[2])
                                {
                                    checkedPosition[2] = couplingBoundMax_[2] - boundCorr_;
                                }
                            }
                        }

                        const label typeIndex = findIndex(typeNames_, rcvParcType[ifacepts][pts]);

                        insertParcel(checkedPosition, velocity, typeIndex);
                        parcelAdded = true;
                        inserted++;
                    }
                }
            }
        }

        if(inserted != 0)
        {
            std::cout << "    Coupling parcels inserted       = " << inserted << std::endl;
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
{}

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

void dsmcMdCoupling::forget(bool forget)
{
	forAll(recvInterfaces_, iface)
	{
		recvInterfaces_[iface]->forget(currIteration_, forget);
	}
}

void dsmcMdCoupling::forget(label iteration, bool forget)
{
	forAll(recvInterfaces_, iface)
	{
		recvInterfaces_[iface]->forget(iteration, forget);
	}
}

void dsmcMdCoupling::forget(label iteration, label interface, bool forget)
{
	recvInterfaces_[interface]->forget(iteration, forget);
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
    /*
    else
    {
        std::cout << "dsmcMdCoupling::insertParcel(): Parcel insertion attempted outside of mesh, parcel not inserted" << std::endl;
    }
    */
}

} // End namespace Foam

// ************************************************************************* //
