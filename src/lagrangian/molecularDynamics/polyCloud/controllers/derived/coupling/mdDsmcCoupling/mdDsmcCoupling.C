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
    sendingRegion_(false),
    receivingRegion_(false),
    sendingBound_(false),
    receivingBound_(false),
    idList(molCloud_.cP().molIds()),
    rU_(molCloud_.redUnits()),
    couplingBounds_(false),
    couplingRegion_(false),
    couplingRegionMin_(vector::zero),
    couplingRegionMax_(vector::zero),
    couplingBoundMin_(vector::zero),
    couplingBoundMax_(vector::zero),
    couplingBoundNorm_(vector::zero),
    couplingBoundZeroThick_(vector(-1, -1, -1)),
    overlapEnergyLimit_(molCloud_.pot().potentialEnergyLimit()),
    overlapIterations_(100),
    prevMolCount_(0),
    currIteration_(0),
    boundCorr_(0),
    meshMin_(VGREAT, VGREAT, VGREAT),
    meshMax_(-VSMALL, -VSMALL, -VSMALL),
    initTemperature_(-VSMALL),
    initKe_(-VSMALL),
    initTemperatureDSMC_(-VSMALL),
    initKeDSMC_(-VSMALL),
    nparcsRcv_(0),
    initScaling_(false),
    scaleType_(0)
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
        sendingRegion_ = true;
        sendingBound_ = true;
    }

    if(recvInterfaces_.size() != 0)
    {
        receivingRegion_ = true;
        receivingBound_ = true;
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
        couplingRegionMin_ /= rU_.refLength();
        regionMinFound = true;
	}

    if (propsDict_.found("couplingRegionMax"))
    {
        couplingRegionMax_ = propsDict_.lookup("couplingRegionMax");
        couplingRegionMax_ /= rU_.refLength();
        regionMaxFound = true;
    }

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

        //- There is an overlap between the coupling boundary and the local mesh so should have found at least 1 intersecting cell
        if(!overlap)
        {
            sendingRegion_ = false;
            receivingRegion_ = false;
        }
    }

    bool boundMinFound = false;
    bool boundMaxFound = false;
    bool boundNormFound = false;

    if (propsDict_.found("couplingBoundMin"))
    {
        couplingBoundMin_ = propsDict_.lookup("couplingBoundMin");
        couplingBoundMin_ /= rU_.refLength();
        boundMinFound = true;
    }

    if (propsDict_.found("couplingBoundMax"))
    {
        couplingBoundMax_ = propsDict_.lookup("couplingBoundMax");
        couplingBoundMax_ /= rU_.refLength();
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
            FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
                         << "A coupling boundary must have a zero thickness in at least one direction"
                         << exit(FatalError);
        }

        boundNormFound = true;
    }

    if((boundMinFound && !boundMaxFound && !boundNormFound) ||
       (boundMaxFound && !boundMinFound && !boundNormFound) ||
       (boundNormFound && !boundMinFound && !boundMaxFound))
    {
        FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
                     << "Cannot find couplingBoundMin, couplingBoundMax and couplingBoundNorm together"
                     << exit(FatalError);
    }
    else
    {
        couplingBounds_ = true;

        const cellList& cells = mesh_.cells();

        vector localMeshExtents = meshMax_ - meshMin_;

        //- Find the whole mesh extents
        vector meshExtents = mesh_.bounds().max() - mesh_.bounds().min();

        // Boundary correction value (0.001% extents) calculated against whole mesh extents for consistency at different parallelisation levels
        vector boundCorr = meshExtents * (1e-3 / 100.0);

        // Pick largest correction value as global in each direction
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

        std::cout << "Boundary correction value: " << boundCorr_ << std::endl;

        point cellMin;
        point cellMax;

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
                                << "Coupling boundary zero thickness in z direction but normal value zero"
                                << exit(FatalError);
               }
            }
        }

        if(intersectingCells_.size() > 0)
        {
            std::cout << "Coupling boundary intersecting cell count: " << intersectingCells_.size() << std::endl;
        }
        else
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
                FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
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

    if (propsDict_.found("potentialEnergyInsertLimit"))
    {
        overlapEnergyLimit_ = readScalar(propsDict_.lookup("potentialEnergyInsertLimit"));
        overlapEnergyLimit_ /= rU_.refEnergy();
    }

    if (propsDict_.found("insertIterations"))
    {
        overlapIterations_ = readLabel(propsDict_.lookup("insertIterations"));
    }

    if (propsDict_.found("output"))
    {
        output_ = Switch(propsDict_.lookup("output"));
    }

    if (propsDict_.found("initialVelScaling"))
    {
        initScaling_ = Switch(propsDict_.lookup("initialVelScaling"));
    }

    if(initScaling_)
    {
        word scalingType
        (
            propsDict_.lookup("scalingType")
        );

        if(scalingType == "temperature" || scalingType == "Temperature")
        {
            scaleType_ = 0;
        }

        if(scalingType == "linearKE" || scalingType == "LinearKE")
        {
            scaleType_ = 1;
        }
    }

#ifdef USE_MUI
    //Initialise exact time sampler for MUI
    chrono_sampler = new mui::chrono_sampler_exact3d();
    rcvPoints_.resize(recvInterfaces_.size());
#endif

    rcvMolType_.resize(recvInterfaces_.size());
    rcvMolId_.resize(recvInterfaces_.size());
    rcvVelX_.resize(recvInterfaces_.size());
    rcvVelY_.resize(recvInterfaces_.size());
    rcvVelZ_.resize(recvInterfaces_.size());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mdDsmcCoupling::~mdDsmcCoupling()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool mdDsmcCoupling::initialConfiguration(label stage)
{
    bool returnVal = false;

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
        if(!molCloud_.cloudVelocityScaled())
        {
            //Calculate initial temperature of whole cloud
            initTemperature_ = calcTemperature();
            std::cout << "Initial MD temperature: " << initTemperature_ << std::endl;

            //Calculate initial KE of whole cloud
            initKe_ = calcAvgLinearKe();
            std::cout << "Initial MD average linear KE per molecule: " << initKe_ << std::endl;

            returnVal = receiveCoupledRegion(true); // Receive ghost molecules in coupled region(s) at time = startTime and commit time=1 to release other side

            //Distribute received temperature from DSMC side to all MPI ranks
            if (Pstream::parRun())
            {
                reduce(initTemperatureDSMC_, maxOp<scalar>());
                reduce(initKeDSMC_, maxOp<scalar>());
            }

            std::cout << "Initial DSMC temperature: " << initTemperatureDSMC_ << std::endl;
            std::cout << "Initial DSMC average linear KE per parcel: " << initKeDSMC_ << std::endl;

            if(initScaling_)
            {
                scalar scaleValue = 0;

                if(scaleType_ == 0) //Temperature
                {
                    if (initTemperature_ > 0)
                    {
                        scaleValue = sqrt(initTemperatureDSMC_ / initTemperature_);
                    }
                }
                else if(scaleType_ == 1)
                {
                    if (initKe_ > 0)
                    {
                        scaleValue = sqrt(initKeDSMC_ / initKe_);
                    }
                }

                //Scale molecule velocity field
                scaleVelocity(scaleValue);

                //Calculate new temperature of whole cloud
                scalar newInitTemperature = calcTemperature();
                std::cout << "Scaled MD temperature: " << newInitTemperature << std::endl;

                //Calculate new KE of whole cloud
                scalar newLinearKE = calcAvgLinearKe();
                std::cout << "Scaled MD average linear KE per molecule: " << newLinearKE << std::endl;
            }
        }
    }

    return returnVal;
}

bool mdDsmcCoupling::controlAfterMove(label stage)
{
	bool returnVal = false;

    if(stage == 1) // Find, delete and send molecules that have passed a coupling boundary
	{
        currIteration_++; // Increment the current iteration

        if(sendingBound_)
        {
            returnVal = findCoupledMolecules(); // Find, collate and delete any molecules that have passed a coupling boundary
        }

        label nmolsSent = sendCoupledMolecules(); // Send any molecules deleted by coupling boundary (non-blocking)

        if(nmolsSent > 0)
        {
            if (Pstream::parRun())
            {
                if(Pstream::master())
                {
                    std::cout << "Coupling boundary molecules pushed: " << nmolsSent << std::endl;
                }
                else
                {
                    std::cout << "[" << time_.time().value() << "s] Coupling boundary molecules pushed: " << nmolsSent << std::endl;
                }
            }
            else
            {
                std::cout << "Coupling boundary molecules pushed: " << nmolsSent << std::endl;
            }
        }
	}
	else if (stage == 2) // Receive any molecules coupling boundary (blocking)
	{
	    nparcsRcv_ = receiveCoupledParcels();
	}
	else if (stage == 3) // Receive ghost molecules in coupled region (blocking)
	{
	    returnVal = receiveCoupledRegion(false);
    }
	else if (stage == 4) // Insert any coupling boundary molecules from stage 2
	{
	    if(receivingBound_)
	    {
	        cplMoleculeInsert nmolsInserted = insertCoupledMolecules();

	        if(nmolsInserted.nmolsInserted > 0)
	        {
	            returnVal = true;

	            if (Pstream::parRun())
                {
                    if(Pstream::master())
                    {
                        std::cout << "Coupling boundary parcels inserted: " << nmolsInserted.nmolsInserted << std::endl;

                        if(nmolsInserted.nIts > 0)
                        {
                            std::cout << "High energy overlaps detected during insertion, total iterations to resolve: " << nmolsInserted.nIts << std::endl;
                        }

                        if(nparcsRcv_ != nmolsInserted.nmolsInserted)
                        {
                            std::cout << "Warning: Number of boundary parcels inserted (" << nmolsInserted.nmolsInserted << ") not equal to number received (" << nparcsRcv_ << ")" << std::endl;
                        }
                    }
                    else
                    {
                        std::cout << "[" << time_.time().value() << "s] Coupling boundary parcels inserted: " << nmolsInserted.nmolsInserted << std::endl;

                        if(nmolsInserted.nIts > 0)
                        {
                            std::cout << "[" << time_.time().value() << "s] High energy overlaps detected during insertion, total iterations to resolve: " << nmolsInserted.nIts << std::endl;
                        }

                        if(nparcsRcv_ != nmolsInserted.nmolsInserted)
                        {
                            std::cout << "[" << time_.time().value() << "s] Warning: Number of boundary parcels inserted (" << nmolsInserted.nmolsInserted << ") not equal to number received (" << nparcsRcv_ << ")" << std::endl;
                        }
                    }
                }
                else
                {
                    std::cout << "Coupling boundary parcels inserted: " << nmolsInserted.nmolsInserted << std::endl;

                    if(nmolsInserted.nIts > 0)
                    {
                        std::cout << "High energy overlaps detected during insertion, total iterations to resolve: " << nmolsInserted.nIts << std::endl;
                    }

                    if(nparcsRcv_ != nmolsInserted.nmolsInserted)
                    {
                        std::cout << "Warning: Number of boundary parcels inserted (" << nmolsInserted.nmolsInserted << ") not equal to number received (" << nparcsRcv_ << ")" << std::endl;
                    }
                }
	        }
	    }
	}

    return returnVal;
}

void mdDsmcCoupling::controlAfterForces()
{
    sendCoupledRegionAcc(); // Send the accelerations for ghost molecules in the coupling region(s) (non-blocking)
}

bool mdDsmcCoupling::receiveCoupledRegion(bool init)
{
    bool molChanged = false;
#ifdef USE_MUI
	if(receivingRegion_)
    {
	    const constantMoleculeProperties& cP = molCloud_.cP();
	    const scalar trackTime = mesh_.time().deltaT().value();
	    label molCount = 0;
	    moleculeInsert newMol;

        // Iterate through all receiving interfaces for this controller and extract a points list
        forAll(recvInterfaces_, iface)
        {
            rcvPoints_[iface].clear();
            rcvMolType_[iface].clear();
            rcvMolId_[iface].clear();
            rcvVelX_[iface].clear();
            rcvVelY_[iface].clear();
            rcvVelZ_[iface].clear();

            //- Extract a list of all molecule locations received from other solver through this interface
            rcvPoints_[iface] = recvInterfaces_[iface]->fetch_points<std::string>("type_region", currIteration_, *chrono_sampler);

            if(rcvPoints_[iface].size() > 0)
            {
                //- Extract a list of all molecule change status values received from other solver through this interface
                rcvMolType_[iface] = recvInterfaces_[iface]->fetch_values<std::string>("type_region", currIteration_, *chrono_sampler);

                //- Extract a list of all molecule Id's received from other solver through this interface
                rcvMolId_[iface] = recvInterfaces_[iface]->fetch_values<label>("id_region", currIteration_, *chrono_sampler);

                //- Extract a list of all molecule velocities received from other solver through this interface
                rcvVelX_[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_x_region", currIteration_, *chrono_sampler);
                rcvVelY_[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_y_region", currIteration_, *chrono_sampler);
                rcvVelZ_[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_z_region", currIteration_, *chrono_sampler);
            }
        }

        //- Go through received values and find any that are not of the type set to be received
        forAll(rcvMolType_, ifacepts)
        {
            if(rcvMolType_[ifacepts].size() > 0)
            {
                std::vector<std::string>::iterator rcvMolTypeIt;
                std::vector<mui::point3d>::iterator rcvPointsIt = rcvPoints_[ifacepts].begin();
                std::vector<label>::iterator rcvMolIdIt = rcvMolId_[ifacepts].begin();
                std::vector<scalar>::iterator rcvVelXIt = rcvVelX_[ifacepts].begin();
                std::vector<scalar>::iterator rcvVelYIt = rcvVelY_[ifacepts].begin();
                std::vector<scalar>::iterator rcvVelZIt = rcvVelZ_[ifacepts].begin();

                for (rcvMolTypeIt = rcvMolType_[ifacepts].begin(); rcvMolTypeIt != rcvMolType_[ifacepts].end(); rcvMolTypeIt++) {
                    const label molId = findIndex(molNames_, *rcvMolTypeIt);

                    if(molId == -1) //- molId not found in local list as one to receive so store it as one to remove from lists
                    {
                        rcvMolType_[ifacepts].erase(rcvMolTypeIt--);
                        rcvPoints_[ifacepts].erase(rcvPointsIt--);
                        rcvMolId_[ifacepts].erase(rcvMolIdIt--);
                        rcvVelX_[ifacepts].erase(rcvVelXIt--);
                        rcvVelY_[ifacepts].erase(rcvVelYIt--);
                        rcvVelZ_[ifacepts].erase(rcvVelZIt--);
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
        forAll(rcvPoints_, ifacepts)
        {
            if(rcvPoints_[ifacepts].size() > 0)
            {
                bool newList = false;
                bool newSize = false;

                //- List size has changed, treat as a new list
                if(molId_[ifacepts].size() != rcvMolId_[ifacepts].size())
                {
                    newList = true;

                    if(molId_[ifacepts].size() == 0) //- Local storage is zero size, so resize here instead of later
                    {
                        //- Resize local storage to new received size
                        molId_[ifacepts].setSize(rcvMolId_[ifacepts].size(), -1);
                        molHistory_[ifacepts].setSize(rcvMolId_[ifacepts].size(), NULL);

                        newSize = false;
                    }
                    else //- Local storage not zero size, so resize later
                    {
                        newSize = true;
                    }
                }
                else //- List size the same, ensure list hasn't changed internally
                {
                    forAll(rcvMolId_[ifacepts], rcvMol)
                    {
                        if(rcvMolId_[ifacepts][rcvMol] != molId_[ifacepts][rcvMol]) //- If a molecule ID changes then treat as new list
                        {
                            newList = true;
                            break; //- Break once first change detected
                        }
                    }
                }

                if(newList)
                {
                    DynamicList<label> molsToDelete;
                    List<polyMolecule*> rcvMolHistory(rcvMolId_[ifacepts].size(), NULL);

                    //- Determine which molecules from the last iteration no longer exist in the received list and so need to be deleted
                    forAll(molId_[ifacepts], currMol)
                    {
                        bool molFound = false;

                        forAll(rcvMolId_[ifacepts], rcvMol)
                        {
                            if(molId_[ifacepts][currMol] == rcvMolId_[ifacepts][rcvMol])
                            {
                                molFound = true;
                                break; //- Break inner loop as match found
                            }
                        }

                        if(!molFound)
                        {
                            molsToDelete.append(molId_[ifacepts][currMol]);
                        }
                    }

                    // Delete any molecules that need to be
                    forAll(molId_[ifacepts], currMol)
                    {
                        forAll(molsToDelete, delMol)
                        {
                            if(molsToDelete[delMol] == molId_[ifacepts][currMol])
                            {
                                if(molHistory_[ifacepts][currMol] != NULL)
                                {
                                    molCloud_.deleteParticle(*molHistory_[ifacepts][currMol]);
                                    molHistory_[ifacepts][currMol] = NULL;
                                    molChanged = true;
                                }
                            }
                        }
                    }

                    molsToDelete.clear();

                    //- Determine which molecules in the received list did not exist in the last iteration and so are new
                    forAll(rcvMolId_[ifacepts], rcvMol)
                    {
                        forAll(molId_[ifacepts], currMol)
                        {
                            if(rcvMolId_[ifacepts][rcvMol] == molId_[ifacepts][currMol])
                            {
                                rcvMolHistory[rcvMol] = molHistory_[ifacepts][currMol];
                                break; //- Break inner loop as match found
                            }
                        }
                    }

                    //- Resize local storage if new list size required it
                    if(newSize)
                    {
                        //- Resize local storage to new received size
                        molId_[ifacepts].setSize(rcvMolId_[ifacepts].size(), -1);
                        molHistory_[ifacepts].setSize(rcvMolId_[ifacepts].size(), NULL);
                    }

                    //- Update the local molecule ID store to the received version
                    forAll(rcvMolId_[ifacepts], mol)
                    {
                        molId_[ifacepts][mol] = rcvMolId_[ifacepts][mol];
                        molHistory_[ifacepts][mol] = rcvMolHistory[mol];
                    }

                    rcvMolHistory.clear();
                }

                for (size_t pts = 0; pts < rcvPoints_[ifacepts].size(); pts++)
                {
                    point checkedPosition((rcvPoints_[ifacepts][pts][0] / rU_.refLength()), (rcvPoints_[ifacepts][pts][1] / rU_.refLength()), (rcvPoints_[ifacepts][pts][2] / rU_.refLength()));

                    label cell = mesh_.findCell(checkedPosition);

                    // Attempt insertion/update only if received molecule falls within local mesh extents
                    if(cell != -1)
                    {
                        const label molId = findIndex(molNames_, rcvMolType_[ifacepts][pts]);

                        vector velocity;
                        velocity[0] = rcvVelX_[ifacepts][pts] / rU_.refVelocity();
                        velocity[1] = rcvVelY_[ifacepts][pts] / rU_.refVelocity();
                        velocity[2] = rcvVelZ_[ifacepts][pts] / rU_.refVelocity();
/*
                        if(couplingRegion_)
                        {
                            bool trunc = false;

                            if(checkedPosition[0] <= couplingRegionMin_[0])
                            {
                                checkedPosition[0] = couplingRegionMin_[0] + boundCorr_;
                                trunc = true;
                            }

                            if(checkedPosition[0] >= couplingRegionMax_[0])
                            {
                                checkedPosition[0] = couplingRegionMax_[0] - boundCorr_;
                                trunc = true;
                            }

                            if(checkedPosition[1] <= couplingRegionMin_[1])
                            {
                                checkedPosition[1] = couplingRegionMin_[1] + boundCorr_;
                                trunc = true;
                            }

                            if(checkedPosition[1] >= couplingRegionMax_[1])
                            {
                                checkedPosition[1] = couplingRegionMax_[1] - boundCorr_;
                                trunc = true;
                            }

                            if(checkedPosition[2] <= couplingRegionMin_[2])
                            {
                                checkedPosition[2] = couplingRegionMin_[2] + boundCorr_;
                                trunc = true;
                            }

                            if(checkedPosition[2] >= couplingRegionMax_[2])
                            {
                                checkedPosition[2] = couplingRegionMax_[2] - boundCorr_;
                                trunc = true;
                            }

                            if(trunc)
                            {
                                std::cout << "receiveCoupledRegion particle boundCorr_ applied" << std::endl;
                            }
                        }
*/
                        if(molHistory_[ifacepts][pts] == NULL) //- This molecule is new so insert it
                        {
                            newMol = insertMolecule(checkedPosition, molIds_[molId], true, velocity);
                            molHistory_[ifacepts][pts] = newMol.mol;
                            molChanged = true;
                        }
                        else //- This molecule already exists, so just update its values
                        {
                            molHistory_[ifacepts][pts]->position()[0] = checkedPosition[0];
                            molHistory_[ifacepts][pts]->position()[1] = checkedPosition[1];
                            molHistory_[ifacepts][pts]->position()[2] = checkedPosition[2];
                            molHistory_[ifacepts][pts]->specialPosition()[0] = checkedPosition[0];
                            molHistory_[ifacepts][pts]->specialPosition()[1] = checkedPosition[1];
                            molHistory_[ifacepts][pts]->specialPosition()[2] = checkedPosition[2];

                            molHistory_[ifacepts][pts]->v()[0] = velocity[0];
                            molHistory_[ifacepts][pts]->v()[1] = velocity[1];
                            molHistory_[ifacepts][pts]->v()[2] = velocity[2];

                            molHistory_[ifacepts][pts]->updateAfterMove(cP, trackTime);
                        }

                        molCount++;
                    }
                }
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
    }

	//- This is the initialisation call so commit at iteration = 1 to allow other side to continue
    if(init)
    {
        if(receivingRegion_)
        {
            forAll(recvInterfaces_, iface)
            {
                //Fetch initial DSMC system temperature
                initTemperatureDSMC_ = recvInterfaces_[iface]->fetch<scalar>("init_temp");
                //Fetch initial DSMC system linear KE
                initKeDSMC_ = recvInterfaces_[iface]->fetch<scalar>("init_ke");
                //Commit at t=1 to allow other side to continue
                recvInterfaces_[iface]->commit(currIteration_);
            }
        }
    }
#endif
    return molChanged;
}

void mdDsmcCoupling::sendCoupledRegionAcc()
{
#ifdef USE_MUI
    if(sendingRegion_)
    {
        polyMolecule* molecule = NULL;

        // Iterate through all sending interfaces for this controller
        forAll(sendInterfaces_, iface)
        {
            if(molHistory_[iface].size() > 0)
            {
                forAll(molHistory_[iface], mol) // Iterate through molecules received by this controller
                {
                    molecule = molHistory_[iface][mol];

                    if(molecule != NULL)
                    {
                        //- Determine whether parcel is of a type set to send by this controller
                        const label typeIndex = findIndex(molNames_, molCloud_.cP().molIds()[molecule->id()]);

                        if(typeIndex != -1)
                        {
                            vector siteForcesAccum(vector::zero);

                            forAll(molecule->siteForces(), s)
                            {
                                siteForcesAccum[0] += molecule->siteForces()[s][0];
                                siteForcesAccum[1] += molecule->siteForces()[s][1];
                                siteForcesAccum[2] += molecule->siteForces()[s][2];
                            }

                            if(siteForcesAccum[0] != 0 || siteForcesAccum[1] != 0 || siteForcesAccum[2] != 0)
                            {
                                const scalar& mass = molCloud_.cP().mass(molecule->id());

                                // Get the molecule centre
                                mui::point3d molCentre;
                                molCentre[0] = molecule->position()[0] * rU_.refLength();
                                molCentre[1] = molecule->position()[1] * rU_.refLength();
                                molCentre[2] = molecule->position()[2] * rU_.refLength();

                                // Push molecule type
                                sendInterfaces_[iface]->push("type_region", molCentre, static_cast<std::string>(molNames_[typeIndex]));

                                // Push molecule ID from receive history
                                sendInterfaces_[iface]->push("id_region", molCentre, static_cast<label>(molId_[iface][mol]));

                                // Push the molecule acceleration to the interface
                                vector acc;
                                acc[0] = ((siteForcesAccum[0] * rU_.refForce()) / (mass * rU_.refMass()));
                                acc[1] = ((siteForcesAccum[1] * rU_.refForce()) / (mass * rU_.refMass()));
                                acc[2] = ((siteForcesAccum[2] * rU_.refForce()) / (mass * rU_.refMass()));

                                sendInterfaces_[iface]->push("acc_x_region", molCentre, acc[0]);
                                sendInterfaces_[iface]->push("acc_y_region", molCentre, acc[1]);
                                sendInterfaces_[iface]->push("acc_z_region", molCentre, acc[2]);
                            }
                        }
                    }
                }
            }

            // Commit (transmit) values to the MUI interface
            sendInterfaces_[iface]->commit(currIteration_);
        }
    }
#endif
}

bool mdDsmcCoupling::findCoupledMolecules()
{
    bool molsDeleted = false;

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

                    if(couplingBoundZeroThick_[0] == 1)
                    {
                        if(couplingBoundNorm_[0] > 0) //- Boundary is positive facing in the x
                        {
                            if(molsInCell[molecule]->position()[0] <= couplingBoundMin_[0])
                            {
                                removeMolecule = true;
                            }
                        }
                        else if (couplingBoundNorm_[0] < 0) //- Boundary is negative facing in the x
                        {
                            if(molsInCell[molecule]->position()[0] >= couplingBoundMax_[0])
                            {
                                removeMolecule = true;
                            }
                        }
                    }

                    if(couplingBoundZeroThick_[1] == 1)
                    {
                        if(couplingBoundNorm_[1] > 0) //- Boundary is positive facing in the y
                        {
                            if(molsInCell[molecule]->position()[1] <= couplingBoundMin_[1])
                            {
                                removeMolecule = true;
                            }
                        }
                        else if (couplingBoundNorm_[1] < 0) //- Boundary is negative facing in the x
                        {
                            if(molsInCell[molecule]->position()[1] >= couplingBoundMax_[1])
                            {
                                removeMolecule = true;
                            }
                        }
                    }

                    if(couplingBoundZeroThick_[2] == 1)
                    {
                        if(couplingBoundNorm_[2] > 0) //- Boundary is positive facing in the z
                        {
                            if(molsInCell[molecule]->position()[2] <= couplingBoundMin_[2])
                            {
                                removeMolecule = true;
                            }
                        }
                        else if (couplingBoundNorm_[2] < 0) //- Boundary is negative facing in the x
                        {
                            if(molsInCell[molecule]->position()[2] >= couplingBoundMax_[2])
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

                        if(couplingBounds_)
                        {
                            if(couplingBoundZeroThick_[0] == 1) //- Boundary has zero thickness in the x
                            {
                                newMolToSend.position[0] = couplingBoundMin_[0];
                            }
                            else
                            {
                                if(newMolToSend.position[0] <= couplingBoundMin_[0])
                                {
                                    newMolToSend.position[0] = couplingBoundMin_[0];
                                }

                                if(newMolToSend.position[0] >= couplingBoundMax_[0])
                                {
                                    newMolToSend.position[0] = couplingBoundMax_[0];
                                }
                            }

                            if(couplingBoundZeroThick_[1] == 1) //- Boundary has zero thickness in the y
                            {
                                newMolToSend.position[1] = couplingBoundMin_[1];
                            }
                            else
                            {
                                if(newMolToSend.position[1] <= couplingBoundMin_[1])
                                {
                                    newMolToSend.position[1] = couplingBoundMin_[1];
                                }

                                if(newMolToSend.position[1] >= couplingBoundMax_[1])
                                {
                                    newMolToSend.position[1] = couplingBoundMax_[1];
                                }
                            }

                            if(couplingBoundZeroThick_[2] == 1) //- Boundary has zero thickness in the z
                            {
                                newMolToSend.position[2] = couplingBoundMin_[2];
                            }
                            else
                            {
                                if(newMolToSend.position[2] <= couplingBoundMin_[2])
                                {
                                    newMolToSend.position[2] = couplingBoundMin_[2];
                                }

                                if(newMolToSend.position[2] >= couplingBoundMax_[2])
                                {
                                    newMolToSend.position[2] = couplingBoundMax_[2];
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
        molCloud_.deleteParticle(*molsToRemove[molecule]);
        molsDeleted = true;
    }

    molsToRemove.clear();

    return molsDeleted;
}

label mdDsmcCoupling::sendCoupledMolecules()
{
    label nmolsSent = 0;
#ifdef USE_MUI
    if(sendingBound_)
    {
        forAll(sendInterfaces_, iface)
        {
            if(molsToSend_.size() != 0)
            {
                forAll(molsToSend_, mols)
                {
                    // Get the molecule centre
                    mui::point3d molCentre;
                    molCentre[0] = molsToSend_[mols].position[0] * rU_.refLength();
                    molCentre[1] = molsToSend_[mols].position[1] * rU_.refLength();
                    molCentre[2] = molsToSend_[mols].position[2] * rU_.refLength();

                    // Push molecule type
                    sendInterfaces_[iface]->push("type_bound", molCentre, static_cast<std::string>(molsToSend_[mols].molType));

                    // Push molecule velocity
                    sendInterfaces_[iface]->push("vel_x_bound", molCentre, molsToSend_[mols].velocity[0] * rU_.refVelocity());
                    sendInterfaces_[iface]->push("vel_y_bound", molCentre, molsToSend_[mols].velocity[1] * rU_.refVelocity());
                    sendInterfaces_[iface]->push("vel_z_bound", molCentre, molsToSend_[mols].velocity[2] * rU_.refVelocity());

                    nmolsSent++;
                }
            }

            // Commit (transmit) values to the coupling interface
            sendInterfaces_[iface]->commit(currIteration_);
        }
    }
#endif
    //- Clear the sent molecules
    molsToSend_.clear();

    return nmolsSent;
}

label mdDsmcCoupling::receiveCoupledParcels()
{
    label nparcsRcv = 0;
#ifdef USE_MUI
    if(receivingBound_)
    {
        // Iterate through all receiving interfaces for this controller and extract a points list for each molecule type handled
        forAll(recvInterfaces_, iface)
        {
            rcvPoints_[iface].clear();
            rcvMolType_[iface].clear();
            rcvVelX_[iface].clear();
            rcvVelY_[iface].clear();
            rcvVelZ_[iface].clear();

            //- Extract a list of all molecule locations
            rcvPoints_[iface] = recvInterfaces_[iface]->fetch_points<std::string>("type_bound", currIteration_, *chrono_sampler);

            if(rcvPoints_[iface].size() > 0)
            {
                //- Extract a list of all molecule types
                rcvMolType_[iface] = recvInterfaces_[iface]->fetch_values<std::string>("type_bound", currIteration_, *chrono_sampler);

                //- Extract a list of all molecule velocities received from other solver through this interface
                rcvVelX_[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_x_bound", currIteration_, *chrono_sampler);
                rcvVelY_[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_y_bound", currIteration_, *chrono_sampler);
                rcvVelZ_[iface] = recvInterfaces_[iface]->fetch_values<scalar>("vel_z_bound", currIteration_, *chrono_sampler);
            }
        }

        //- Go through received values and find any that are not of the type set to be received
        forAll(rcvMolType_, ifacepts)
        {
            if(rcvMolType_[ifacepts].size() > 0)
            {
                std::vector<std::string>::iterator rcvMolTypeIt;
                std::vector<mui::point3d>::iterator rcvPointsIt = rcvPoints_[ifacepts].begin();
                std::vector<scalar>::iterator rcvVelXIt = rcvVelX_[ifacepts].begin();
                std::vector<scalar>::iterator rcvVelYIt = rcvVelY_[ifacepts].begin();
                std::vector<scalar>::iterator rcvVelZIt = rcvVelZ_[ifacepts].begin();

                for (rcvMolTypeIt = rcvMolType_[ifacepts].begin(); rcvMolTypeIt != rcvMolType_[ifacepts].end(); rcvMolTypeIt++) {
                    const label molId = findIndex(molNames_, *rcvMolTypeIt);

                    if(molId == -1) //- molId not found in local list as one to receive so store it as one to remove from lists
                    {
                        rcvMolType_[ifacepts].erase(rcvMolTypeIt--);
                        rcvPoints_[ifacepts].erase(rcvPointsIt--);
                        rcvVelX_[ifacepts].erase(rcvVelXIt--);
                        rcvVelY_[ifacepts].erase(rcvVelYIt--);
                        rcvVelZ_[ifacepts].erase(rcvVelZIt--);
                    }

                    rcvPointsIt++;
                    rcvVelXIt++;
                    rcvVelYIt++;
                    rcvVelZIt++;
                }
            }
        }

        // Iterate through points received
        forAll(rcvPoints_, ifacepts)
        {
            if(rcvPoints_[ifacepts].size() > 0)
            {
                for (size_t pts = 0; pts < rcvPoints_[ifacepts].size(); pts++)
                {
                    point checkedPosition(rcvPoints_[ifacepts][pts][0] / rU_.refLength(), rcvPoints_[ifacepts][pts][1] / rU_.refLength(), rcvPoints_[ifacepts][pts][2] / rU_.refLength());

                    label cell = mesh_.findCell(checkedPosition);

                    // Attempt insertion/update only if received molecule falls within local mesh extents
                    if(cell != -1)
                    {
                        vector velocity;
                        velocity[0] = rcvVelX_[ifacepts][pts] / rU_.refVelocity();
                        velocity[1] = rcvVelY_[ifacepts][pts] / rU_.refVelocity();
                        velocity[2] = rcvVelZ_[ifacepts][pts] / rU_.refVelocity();
/*
                        if(couplingBounds_)
                        {
                            bool trunc = false;

                            if(couplingBoundZeroThick_[0] == 1) //- Boundary has zero thickness in the x
                            {
                                if(checkedPosition[0] <= couplingBoundMin_[0])
                                {
                                    checkedPosition[0] = couplingBoundMin_[0] + (couplingBoundNorm_[0] * boundCorr_);
                                    trunc = true;
                                }
                            }
                            else
                            {
                                if(checkedPosition[0] <= couplingBoundMin_[0])
                                {
                                    checkedPosition[0] = couplingBoundMin_[0] + boundCorr_;
                                    trunc = true;
                                }

                                if(checkedPosition[0] >= couplingBoundMax_[0])
                                {
                                    checkedPosition[0] = couplingBoundMax_[0] - boundCorr_;
                                    trunc = true;
                                }
                            }

                            if(couplingBoundZeroThick_[1] == 1) //- Boundary has zero thickness in the y
                            {
                                if(checkedPosition[1] <= couplingBoundMin_[1])
                                {
                                    checkedPosition[1] = couplingBoundMin_[1] + (couplingBoundNorm_[1] * boundCorr_);
                                    trunc = true;
                                }
                            }
                            else
                            {
                                if(checkedPosition[1] <= couplingBoundMin_[1])
                                {
                                    checkedPosition[1] = couplingBoundMin_[1] + boundCorr_;
                                    trunc = true;
                                }

                                if(checkedPosition[1] >= couplingBoundMax_[1])
                                {
                                    checkedPosition[1] = couplingBoundMax_[1] - boundCorr_;
                                    trunc = true;
                                }
                            }

                            if(couplingBoundZeroThick_[2] == 1) //- Boundary has zero thickness in the z
                            {
                                if(checkedPosition[2] <= couplingBoundMin_[2])
                                {
                                    checkedPosition[2] = couplingBoundMin_[2] + (couplingBoundNorm_[2] * boundCorr_);
                                    trunc = true;
                                }
                            }
                            else
                            {
                                if(checkedPosition[2] <= couplingBoundMin_[2])
                                {
                                    checkedPosition[2] = couplingBoundMin_[2] + boundCorr_;
                                    trunc = true;
                                }

                                if(checkedPosition[2] >= couplingBoundMax_[2])
                                {
                                    checkedPosition[2] = couplingBoundMax_[2] - boundCorr_;
                                    trunc = true;
                                }
                            }

                            if(trunc)
                            {
                                std::cout << "receiveCoupledParcels particle boundCorr_ applied" << std::endl;
                            }
                        }
                        */

                        coupledMolecule newMol;
                        newMol.molType = rcvMolType_[ifacepts][pts];
                        newMol.position = checkedPosition;
                        newMol.velocity = velocity;
                        molsReceived_.append(newMol);

                        nparcsRcv++;
                    }
                }
            }
        }
    }
#endif
    return nparcsRcv;
}

mdDsmcCoupling::cplMoleculeInsert mdDsmcCoupling::insertCoupledMolecules()
{
    cplMoleculeInsert newMolInsert;
#ifdef USE_MUI
    if(molsReceived_.size() > 0)
    {
        moleculeInsert newMol;

        forAll(molsReceived_, mol)
        {
            const label molId = findIndex(molNames_, molsReceived_[mol].molType);

            point molPosition = molsReceived_[mol].position;

            newMol = insertMolecule(molPosition, molIds_[molId], false, molsReceived_[mol].velocity);

            if(newMol.mol != NULL)
            {
                newMolInsert.nmolsInserted++;
                newMolInsert.nIts += newMol.addedIterations;
            }
        }
    }
#endif
    molsReceived_.clear();

    return newMolInsert;
}

void mdDsmcCoupling::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
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

void mdDsmcCoupling::forget(bool forget)
{
	forAll(recvInterfaces_, iface)
	{
		recvInterfaces_[iface]->forget(currIteration_, forget);
	}
}

void mdDsmcCoupling::forget(label iteration, bool forget)
{
	forAll(recvInterfaces_, iface)
	{
		recvInterfaces_[iface]->forget(iteration, forget);
	}
}

void mdDsmcCoupling::forget(label iteration, label interface, bool forget)
{
	recvInterfaces_[interface]->forget(iteration, forget);
}

mdDsmcCoupling::moleculeInsert mdDsmcCoupling::insertMolecule
(
    point& position,
    const label& id,
    const bool ghost,
    vector& velocity
)
{
    label cell = mesh_.findCell(position);
    moleculeInsert newInsert;
    
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

        polyMolecule* newMol = molCloud_.createOnlyMolecule
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

        if (ghost) //Ghost molecule insertion so no overlap testing needed, just insert into cloud
        {
            newInsert.mol = newMol;
            molCloud_.insertMolecule(newMol);

            molCloud_.updateNeighbouringRadii(newMol);
            molCloud_.insertMolInCellOccupancy(newMol);
        }
        else //If molecule not a ghost then perform overlap testing
        {
            polyMolecule* overlapMol = NULL;
            overlapMol = checkForOverlaps(newMol, overlapEnergyLimit_);
            label iterCount = 0;

            if(overlapMol != NULL) // Check for any overlaps that will exceed pre-determined energy limit for interactions
            {
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
                // Limit for maximum perturbation distance calculation equal to quarter of total extents of local mesh
                vector maxPerturb = (meshMax_ - meshMin_) * 0.25;

                for (label i=0; i<overlapIterations_; i++)
                {
                    // Delete the created molecule that is causing the overlap
                    if(newMol != NULL)
                    {
                        delete(newMol);
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

                    //- Calculate perturbed position
                    vector perturbation(perturbDistance * overlapNorm[0], perturbDistance * overlapNorm[1], perturbDistance * overlapNorm[2]);

                    //Ensure perturbation doesn't exceed maximum for the local mesh extents
                    if(perturbation[0] > maxPerturb[0])
                    {
                        perturbation[0] = maxPerturb[0];
                    }

                    if(perturbation[1] > maxPerturb[1])
                    {
                        perturbation[1] = maxPerturb[1];
                    }

                    if(perturbation[2] > maxPerturb[2])
                    {
                        perturbation[2] = maxPerturb[2];
                    }

                    //Add perturbation to molecule position
                    position[0] += perturbation[0];
                    position[1] += perturbation[1];
                    position[2] += perturbation[2];

                    //- If there are coupling boundary details then make sure the perturbed position doesn't fall outside of it
                    if(couplingBounds_)
                    {
                        trunkX = false;
                        trunkY = false;
                        trunkZ = false;

                        if(couplingBoundZeroThick_[0] == 1) //- Boundary has zero thickness in the x
                        {
                            if(couplingBoundNorm_[0] > 0)
                            {
                                if(position[0] <= couplingBoundMin_[0])
                                {
                                    position[0] = couplingBoundMin_[0] + boundCorr_;
                                    trunkX = true;
                                }
                            }
                            else if (couplingBoundNorm_[0] < 0)
                            {
                                if(position[0] >= couplingBoundMax_[0])
                                {
                                    position[0] = couplingBoundMax_[0] - boundCorr_;
                                    trunkX = true;
                                }
                            }
                        }
                        else
                        {
                            if(position[0] <= couplingBoundMin_[0])
                            {
                                position[0] = couplingBoundMin_[0] + boundCorr_;
                                trunkX = true;
                            }

                            if(position[0] >= couplingBoundMax_[0])
                            {
                                position[0] = couplingBoundMax_[0] - boundCorr_;
                                trunkX = true;
                            }
                        }

                        if(couplingBoundZeroThick_[1] == 1) //- Boundary has zero thickness in the y
                        {
                            if(couplingBoundNorm_[1] > 0)
                            {
                                if(position[1] <= couplingBoundMin_[1])
                                {
                                    position[1] = couplingBoundMin_[1] + boundCorr_;
                                    trunkY = true;
                                }
                            }
                            else if(couplingBoundNorm_[1] < 0)
                            {
                                if(position[1] >= couplingBoundMax_[1])
                                {
                                    position[1] = couplingBoundMax_[1] - boundCorr_;
                                    trunkY = true;
                                }
                            }
                        }
                        else
                        {
                            if(position[1] <= couplingBoundMin_[1])
                            {
                                position[1] = couplingBoundMin_[1] + boundCorr_;
                                trunkY = true;
                            }

                            if(position[1] >= couplingBoundMax_[1])
                            {
                                position[1] = couplingBoundMax_[1] - boundCorr_;
                                trunkY = true;
                            }
                        }

                        if(couplingBoundZeroThick_[2] == 1) //- Boundary has zero thickness in the z
                        {
                            if(couplingBoundNorm_[2] > 0)
                            {
                                if(position[2] <= couplingBoundMin_[2])
                                {
                                    position[2] = couplingBoundMin_[2] + boundCorr_;
                                    trunkZ = true;
                                }
                            }
                            else if(couplingBoundNorm_[2] < 0)
                            {
                                if(position[2] >= couplingBoundMax_[2])
                                {
                                    position[2] = couplingBoundMax_[2] - boundCorr_;
                                    trunkZ = true;
                                }
                            }
                        }
                        else
                        {
                            if(position[2] <= couplingBoundMin_[2])
                            {
                                position[2] = couplingBoundMin_[2] + boundCorr_;
                                trunkZ = true;
                            }

                            if(position[2] >= couplingBoundMax_[2])
                            {
                                position[2] = couplingBoundMax_[2] - boundCorr_;
                                trunkZ = true;
                            }
                        }

                        //Check against local mesh extents
                        if(position[0] <= meshMin_[0])
                        {
                            position[0] = meshMin_[0];
                            trunkX = true;
                        }

                        if(position[0] >= meshMax_[0])
                        {
                            position[0] = meshMax_[0];
                            trunkX = true;
                        }

                        if(position[1] <= meshMin_[1])
                        {
                            position[1] = meshMin_[1];
                            trunkY = true;
                        }

                        if(position[1] >= meshMax_[1])
                        {
                            position[1] = meshMax_[1];
                            trunkY = true;
                        }

                        if(position[2] <= meshMin_[2])
                        {
                            position[2] = meshMin_[2];
                            trunkZ = true;
                        }

                        if(position[2] >= meshMax_[2])
                        {
                            position[2] = meshMax_[2];
                            trunkZ = true;
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
                        //- Undo the initial perturbation so new randomised normal can be used instead
                        position = origPos;

                        overlapNorm[0] = molCloud_.rndGen().position<scalar>(-1.0, 1.0);
                        overlapNorm[1] = molCloud_.rndGen().position<scalar>(-1.0, 1.0);
                        overlapNorm[2] = molCloud_.rndGen().position<scalar>(-1.0, 1.0);

                        //- Calculate perturbed position
                        vector perturbation(perturbDistance * overlapNorm[0], perturbDistance * overlapNorm[1], perturbDistance * overlapNorm[2]);

                        //Ensure perturbation doesn't exceed maximum for the local mesh extents
                        if(perturbation[0] > maxPerturb[0])
                        {
                            perturbation[0] = maxPerturb[0];
                        }

                        if(perturbation[1] > maxPerturb[1])
                        {
                            perturbation[1] = maxPerturb[1];
                        }

                        if(perturbation[2] > maxPerturb[2])
                        {
                            perturbation[2] = maxPerturb[2];
                        }

                        //Add perturbation to molecule position using new randomised normal
                        position[0] += perturbation[0];
                        position[1] += perturbation[1];
                        position[2] += perturbation[2];

                        //- If there are coupling boundary details then make sure the perturbed position doesn't fall outside of it
                        if(couplingBounds_)
                        {
                            if(couplingBoundZeroThick_[0] == 1) //- Boundary has zero thickness in the x
                            {
                                if(couplingBoundNorm_[0] > 0)
                                {
                                    if(position[0] <= couplingBoundMin_[0])
                                    {
                                        position[0] = couplingBoundMin_[0] + boundCorr_;
                                    }
                                }
                                else if (couplingBoundNorm_[0] < 0)
                                {
                                    if(position[0] >= couplingBoundMax_[0])
                                    {
                                        position[0] = couplingBoundMax_[0] - boundCorr_;
                                    }
                                }
                            }
                            else
                            {
                                if(position[0] <= couplingBoundMin_[0])
                                {
                                    position[0] = couplingBoundMin_[0] + boundCorr_;
                                }

                                if(position[0] >= couplingBoundMax_[0])
                                {
                                    position[0] = couplingBoundMax_[0] - boundCorr_;
                                }
                            }

                            if(couplingBoundZeroThick_[1] == 1) //- Boundary has zero thickness in the y
                            {
                                if(couplingBoundNorm_[1] > 0)
                                {
                                    if(position[1] <= couplingBoundMin_[1])
                                    {
                                        position[1] = couplingBoundMin_[1] + boundCorr_;
                                    }
                                }
                                else if(couplingBoundNorm_[1] < 0)
                                {
                                    if(position[1] >= couplingBoundMax_[1])
                                    {
                                        position[1] = couplingBoundMax_[1] - boundCorr_;
                                    }
                                }
                            }
                            else
                            {
                                if(position[1] <= couplingBoundMin_[1])
                                {
                                    position[1] = couplingBoundMin_[1] + boundCorr_;
                                }

                                if(position[1] >= couplingBoundMax_[1])
                                {
                                    position[1] = couplingBoundMax_[1] - boundCorr_;
                                }
                            }

                            if(couplingBoundZeroThick_[2] == 1) //- Boundary has zero thickness in the z
                            {
                                if(couplingBoundNorm_[2] > 0)
                                {
                                    if(position[2] <= couplingBoundMin_[2])
                                    {
                                        position[2] = couplingBoundMin_[2] + boundCorr_;
                                    }
                                }
                                else if(couplingBoundNorm_[2] < 0)
                                {
                                    if(position[2] >= couplingBoundMax_[2])
                                    {
                                        position[2] = couplingBoundMax_[2] - boundCorr_;
                                    }
                                }
                            }
                            else
                            {
                                if(position[2] <= couplingBoundMin_[2])
                                {
                                    position[2] = couplingBoundMin_[2] + boundCorr_;
                                }

                                if(position[2] >= couplingBoundMax_[2])
                                {
                                    position[2] = couplingBoundMax_[2] - boundCorr_;
                                }
                            }

                            //Check against local mesh extents
                            if(position[0] <= meshMin_[0])
                            {
                                position[0] = meshMin_[0];
                            }

                            if(position[0] >= meshMax_[0])
                            {
                                position[0] = meshMax_[0];
                            }

                            if(position[1] <= meshMin_[1])
                            {
                                position[1] = meshMin_[1];
                            }

                            if(position[1] >= meshMax_[1])
                            {
                                position[1] = meshMax_[1];
                            }

                            if(position[2] <= meshMin_[2])
                            {
                                position[2] = meshMin_[2];
                            }

                            if(position[2] >= meshMax_[2])
                            {
                                position[2] = meshMax_[2];
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
                        newMol = molCloud_.createOnlyMolecule
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
                            if(newMol != NULL)
                            {
                                // Delete the created molecule as it will exceed energy limit according to force-field calculation
                                delete(newMol);
                            }

                            std::cout << "mdDsmcCoupling::insertMolecule(): Molecule insertion attempted outside of mesh whilst finding new location, molecule not inserted" << std::endl;
                            std::cout << "Attempted pos: " << position[0] << "," << position[1] << "," << position[2] << std::endl;
                            return newInsert;
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

            if(overlapMol != NULL) //The iterative process to perturb the molecule away from overlap failed after nIter tries, inform and return null
            {
                std::cout << "mdDsmcCoupling::insertMolecule(): Failed to find new location for molecule, not inserted" << std::endl;

                // Delete offending molecule before returning
                if(newMol != NULL)
                {
                    delete(newMol);
                }

                return newInsert;
            }
            else if (iterCount > 0) //Overlapping molecule(s) found but corrected in a number of iterations
            {
                newInsert.addedIterations = iterCount;
                newInsert.mol = newMol;
                molCloud_.insertMolecule(newMol);

                molCloud_.updateNeighbouringRadii(newMol);
                molCloud_.insertMolInCellOccupancy(newMol);
            }
            else //No overlapping molecule found and no iterations needed
            {
                newInsert.mol = newMol;
                molCloud_.insertMolecule(newMol);

                molCloud_.updateNeighbouringRadii(newMol);
                molCloud_.insertMolInCellOccupancy(newMol);
            }
        }

        //Molecule insertion succeeded, return pointer to new molecule
        return newInsert;
    }
    else //Failed to find the cell in the mesh to insert (shouldn't happen), return null
    {
        return newInsert;
    }
}

polyMolecule* mdDsmcCoupling::checkForOverlaps(polyMolecule* newMol, const scalar& potEnergyLimit)
{
	polyMolecule* molJ = NULL;

	// Real-Real interactions
	const cellInteractions<polyMolecule>& il = molCloud_.il();
	const labelListList& dil = il.dil();
	const List<DynamicList<polyMolecule*> >& cellOccupancy = molCloud_.cellOccupancy();

	{
        forAll(dil, d)
        {
            forAll(dil[d], interactingCells)
            {
                List<polyMolecule*> cellJ =	cellOccupancy[dil[d][interactingCells]];

                forAll(cellJ, cellJMols)
                {
                    molJ = cellJ[cellJMols];

                    if(molCloud_.evaluatePotentialLimit(newMol, molJ, potEnergyLimit))
                    {
                        return molJ;
                    }
                }
            }

            forAll(cellOccupancy[d], cellIOtherMols)
            {
                molJ = cellOccupancy[d][cellIOtherMols];

                if(molCloud_.evaluatePotentialLimit(newMol, molJ, potEnergyLimit))
                {
                    return molJ;
                }
            }
        }
	}

	{
        // Real-Referred interactions
        forAll(il.refCellsParticles(), r)
        {
            forAll(il.refCellsParticles()[r], i)
            {
                molJ = il.refCellsParticles()[r][i];

                if(molCloud_.evaluatePotentialLimit(newMol, molJ, potEnergyLimit))
                {
                    return molJ;
                }
            }
        }
	}

	return NULL;
}

scalar mdDsmcCoupling::calcTemperature()
{
    // - calculate streaming velocity
    scalar mass = 0;
    vector mom = vector::zero;
    scalar kE = 0;
    scalar dof = 0;
    scalar angularKeSum = 0;

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for(mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(!mol().ghost() && findIndex(molIds_, mol().id()) != -1)
        {
            const scalar& massI = molCloud_.cP().mass(mol().id());

            mass += massI;
            mom += massI*mol().v();
        }
    }

    if (Pstream::parRun())
    {
        reduce(mom, sumOp<vector>());
        reduce(mass, sumOp<scalar>());
    }

    vector velocity(vector::zero);

    if(mass > 0)
    {
        velocity = mom/mass;
    }

    for(mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(!mol().ghost() && findIndex(molIds_, mol().id()) != -1)
        {
            const scalar& massI = molCloud_.cP().mass(mol().id());

            kE += 0.5*massI*magSqr(mol().v() - velocity);
            dof += molCloud_.cP().degreesOfFreedom(mol().id());

            const diagTensor& molMoI(molCloud_.cP().momentOfInertia(mol().id()));

            // angular speed
            const vector& molOmega(inv(molMoI) & mol().pi());
            angularKeSum += 0.5*(molOmega & molMoI & molOmega);
        }
    }

    if (Pstream::parRun())
    {
        reduce(kE, sumOp<scalar>());
        reduce(dof, sumOp<scalar>());
        reduce(angularKeSum, sumOp<scalar>());
    }

    scalar tempMeasI = 0;

    if(dof > 0.0)
    {
        const scalar& kB = molCloud_.redUnits().kB();

        tempMeasI = (2.0*(kE+angularKeSum))/(kB*dof);

        return tempMeasI*rU_.refTemp();
    }
    else
    {
       return tempMeasI;
    }
}

scalar mdDsmcCoupling::calcAvgLinearKe()
{
    scalar avgKe = 0;
    label count = 0;

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for(mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(!mol().ghost() && findIndex(molIds_, mol().id()) != -1)
        {
            const scalar& massI = molCloud_.cP().mass(mol().id()) * rU_.refMass();
            avgKe += (0.5 * massI)*(magSqr(mol().v() * rU_.refVelocity()));
            count++;
        }
    }

    if (Pstream::parRun())
    {
        reduce(avgKe, sumOp<scalar>());
        reduce(count, sumOp<label>());
    }

    if(avgKe > 0)
    {
        avgKe /= static_cast<scalar>(count);
    }

    return avgKe;
}

void mdDsmcCoupling::scaleVelocity(scalar scaleValue)
{
    std::cout << "Scaling velocity using (" << scaleValue << ")" << std::endl;

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for(mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(!mol().ghost())
        {
            mol().v() *= scaleValue;
        }
    }

    molCloud_.setVelocityScaled();
}

} // End namespace Foam

// ************************************************************************* //
