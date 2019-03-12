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
    output_(false),
#ifdef USE_MUI
    cellCentres_(),
    sendInterfaces_(),
    recvInterfaces_(),
#endif
    sending_(false),
    receiving_(false)
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
		//Initialise exact time and spatial samplers for MUI with a numerical tolerance of 1e-9, large value needed as dsmcFoamPlus works in non-normalised numerics.
		spatial_sampler = new mui::sampler_exact3d<scalar>(1e-9);
		chrono_sampler = new mui::chrono_sampler_exact3d(1e-9);
#endif
    }

    meshMin_[0] = VGREAT;
	meshMin_[1] = VGREAT;
	meshMin_[2] = VGREAT;
	meshMax_[0] = -VSMALL;
	meshMax_[1] = -VSMALL;
	meshMax_[2] = -VSMALL;

	const pointField& meshPoints = mesh_.points();

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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcMdCoupling::~dsmcMdCoupling()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcMdCoupling::initialConfiguration()
{
    if(sending_)
    {
        sendCoupledRegion(true); // Send ghost parcels in coupled regions at time = startTime
        barrier(static_cast<scalar>(1.0));
    }
}

void dsmcMdCoupling::controlParcelsBeforeCollisions()
{
	if(receiving_)
	{
		if(receiveCoupledMolecules()) // Receive any molecules from MD coupling boundary (blocking)
		{
			cloud_.reBuildCellOccupancy(); // Rebuild cell occupancy if parcels received
		}
	}

	if(sending_)
	{
		sendCoupledParcels(); // Send any parcels deleted by coupling boundary (non-blocking)
	}
}

void dsmcMdCoupling::calculateProperties()
{
    if(sending_)
    {
    	sendCoupledRegion(false); // Send ghost molecules in coupled regions (non-blocking)
    }
}

void dsmcMdCoupling::sendCoupledRegion(bool init)
{
#ifdef USE_MUI
    dsmcParcel* parcel = NULL;
    scalar couplingTime;
    if(init)
    {
    	couplingTime = 1.0;
    }
    else
    {
    	couplingTime = time_.time().value() * oneOverRefTime_;
    }
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
				const label& cellI = controlZone(regionIds()[id])[c];
				const List<dsmcParcel*>& parcelsInCell = cloud_.cellOccupancy()[cellI];

				forAll(parcelsInCell, p) //Iterate through parcels in cell to determine if list size has changed
				{
					parcel = parcelsInCell[p];
					const label typeIndex = findIndex(cloud_.typeIdList(), cloud_.typeIdList()[parcel->typeId()]);

					if(typeIndex != -1)
					{
						currParcelsInCells.append(parcel->origId());
					}
				}
			}
		}

		if(parcelsInCellHistory_[iface].size() != currParcelsInCells.size())
		{
			parcelsInCellHistory_[iface].clear();
			parcelsInCellHistory_[iface].setSize(currParcelsInCells.size());
			listSizeChanged[iface] = true;
		}

		std::cout << "Number of molecules sent to coupled region = " << currParcelsInCells.size() << std::endl;

		forAll(regionIds(), id)
		{
			forAll(controlZone(regionIds()[id]), c)
			{
				const label& cellI = controlZone(regionIds()[id])[c];
				const List<dsmcParcel*>& parcelsInCell = cloud_.cellOccupancy()[cellI];

				forAll(parcelsInCell, p) // Iterate through parcels in cell
				{
					parcel = parcelsInCell[p];

					const label typeIndex = findIndex(cloud_.typeIdList(), cloud_.typeIdList()[parcel->typeId()]);

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

						sendStr.str("");
						sendStr.clear();
						sendStr << typeNames_[typeIndex] << "_" << "parc_changed"; //Send string in format [type]_parc_changed

						// Send flag to say whether this parcel is new or not (to avoid creating mdFoamPlus molecules where they don't have to be)
						sendInterfaces_[iface]->push(sendStr.str(), parcCentre, static_cast<scalar>(parcelChanged));

						// Push the parcel velocity to the interface
						sendInterfaces_[iface]->push("parc_vel_x_region", parcCentre, parcel->U()[0]);
						sendInterfaces_[iface]->push("parc_vel_y_region", parcCentre, parcel->U()[1]);
						sendInterfaces_[iface]->push("parc_vel_z_region", parcCentre, parcel->U()[2]);
					}
				}
			}
		}
	}

	forAll(sendInterfaces_, iface)
    {
		// Commit (transmit) values to the MUI interface
        sendInterfaces_[iface]->commit(couplingTime);
    }
#endif
}

void dsmcMdCoupling::sendCoupledParcels()
{
#ifdef USE_MUI
	const DynamicList<dsmcCloud::coupledParcs>& parcsToSend = cloud_.coupledParcels();
	scalar couplingTime = time_.time().value() * oneOverRefTime_;
	std::stringstream sendStr;

	if(parcsToSend.size() != 0)
	{
		forAll(parcsToSend, parcs)
		{
			const label typeIndex = findIndex(typeNames_, parcsToSend[parcs].parcType);

			if(typeIndex != -1)
			{
				forAll(parcsToSend[parcs].sendingInterfaces, interface)
				{
					const label iface = findIndex(sendInterfaceNames_, parcsToSend[parcs].sendingInterfaces[interface]);

					if(iface != -1)
					{
						// Get the parcel centre
						mui::point3d parcCentre;
						parcCentre[0] = parcsToSend[parcs].parcel->position()[0] * oneOverRefLength_;
						parcCentre[1] = parcsToSend[parcs].parcel->position()[1] * oneOverRefLength_;
						parcCentre[2] = parcsToSend[parcs].parcel->position()[2] * oneOverRefLength_;

						// Push the molecule velocity to the interface
						sendStr.str("");
						sendStr.clear();
						sendStr << parcsToSend[parcs].parcType << "_" << "parc_vel_x_bound";
						sendInterfaces_[iface]->push(sendStr.str(), parcCentre, parcsToSend[parcs].parcel->U()[0]);

						sendStr.str("");
						sendStr.clear();
						sendStr << parcsToSend[parcs].parcType << "_" << "parc_vel_y_bound";
						sendInterfaces_[iface]->push(sendStr.str(), parcCentre, parcsToSend[parcs].parcel->U()[1]);

						sendStr.str("");
						sendStr.clear();
						sendStr << parcsToSend[parcs].parcType << "_" << "parc_vel_z_bound";
						sendInterfaces_[iface]->push(sendStr.str(), parcCentre, parcsToSend[parcs].parcel->U()[2]);

						std::cout << "Coupling boundary parcel pushed at: " << "[" << parcCentre[0] << "," << parcCentre[1] << "," << parcCentre[2] << "]" << std::endl;
					}
				}
			}
		}
	}

	forAll(sendInterfaces_, iface)
	{
		// Commit (transmit) values to the MUI interface
		sendInterfaces_[iface]->commit(couplingTime);
	}
#endif
}

bool dsmcMdCoupling::receiveCoupledMolecules()
{
	bool parcelAdded = false;
#ifdef USE_MUI
	scalar couplingTime = time_.time().value() * oneOverRefTime_;
	std::vector<mui::point3d> rcvPoints;
	std::stringstream rcvStr;
		
	// Iterate through all receiving interfaces for this controller
	forAll(recvInterfaces_, iface)
	{
		forAll(typeNames_, molType)
		{
			rcvStr.str("");
			rcvStr.clear();
			rcvStr << typeNames_[molType] << "_" << "mol_vel_x_bound"; //Receive string in format [type]_mol_vel_x_bound

			//- Extract a list of all molecule locations received from other solver through this interface
			rcvPoints = recvInterfaces_[iface]->fetch_points<scalar>(rcvStr.str(), couplingTime, *chrono_sampler);

			if(rcvPoints.size() > 0)
			{
				const label typeId = findIndex(cloud_.typeIdList(), typeNames_[molType]);

				for (size_t pts = 0; pts < rcvPoints.size(); pts++)
				{
					std::cout << "Coupling boundary molecule received at: [" << rcvPoints[pts][0] << "," << rcvPoints[pts][1] << "," << rcvPoints[pts][2] << "]" << std::endl;

					vector velocity;

					rcvStr.str("");
					rcvStr.clear();
					rcvStr << typeNames_[molType] << "_" << "mol_vel_x_bound";
					velocity[0] = recvInterfaces_[iface]->fetch(rcvStr.str(), rcvPoints[pts], couplingTime, *spatial_sampler, *chrono_sampler);
					rcvStr.str("");
					rcvStr.clear();
					rcvStr << typeNames_[molType] << "_" << "mol_vel_y_bound";
					velocity[1] = recvInterfaces_[iface]->fetch(rcvStr.str(), rcvPoints[pts], couplingTime, *spatial_sampler, *chrono_sampler);
					rcvStr.str("");
					rcvStr.clear();
					rcvStr << typeNames_[molType] << "_" << "mol_vel_z_bound";
					velocity[2] = recvInterfaces_[iface]->fetch(rcvStr.str(), rcvPoints[pts], couplingTime, *spatial_sampler, *chrono_sampler);

					const point position(rcvPoints[pts][0] * refLength_, rcvPoints[pts][1] * refLength_, rcvPoints[pts][2] * refLength_);
					insertParcel(position, velocity, typeId);
					parcelAdded = true;
				}
			}
		}
	}
#endif
	return parcelAdded;
}

void dsmcMdCoupling::output
(
    const fileName& fixedPathName,
    const List<fileName>& timePaths
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
	scalar barrierTime = chrono_sampler->get_lower_bound(time_.time().value() * oneOverRefTime_);

	forAll(sendInterfaces_, iface)
	{
		sendInterfaces_[iface]->barrier(barrierTime);
	}
}

void dsmcMdCoupling::barrier(scalar time)
{
	forAll(sendInterfaces_, iface)
	{
		sendInterfaces_[iface]->barrier(time);
	}
}

void dsmcMdCoupling::barrier(label interface)
{
	scalar barrierTime = chrono_sampler->get_lower_bound(time_.time().value() * oneOverRefTime_);
	sendInterfaces_[interface]->barrier(barrierTime);
}

void dsmcMdCoupling::barrier(scalar time, label interface)
{
	// Wait for the other side to catch up
	sendInterfaces_[interface]->barrier(time);
}

void dsmcMdCoupling::forget()
{
	scalar time = time_.time().value() * oneOverRefTime_;
	forAll(recvInterfaces_, iface)
	{
		recvInterfaces_[iface]->forget(chrono_sampler->get_upper_bound(time), true);
	}
}

void dsmcMdCoupling::forget(scalar time)
{
	forAll(recvInterfaces_, iface)
	{
		recvInterfaces_[iface]->forget(chrono_sampler->get_upper_bound(time), true);
	}
}

void dsmcMdCoupling::forget(scalar time, label interface)
{
	recvInterfaces_[interface]->forget(chrono_sampler->get_upper_bound(time), true);
}

void dsmcMdCoupling::insertParcel
(
    const point& position,
	const vector& U,
	const label& typeId
)
{
	//These need to be properly defined, only correct for sample Argon case
	const scalar rotationalTemperature = 0;
	const scalar vibrationalTemperature = 0;
	const scalar electronicTemperature = 0;

    label cell = -1;
    label tetFace = -1;
    label tetPt = -1;

	mesh_.findCellFacePt
	(
		position,
		cell,
		tetFace,
		tetPt
	);

	if(cell != -1)
	{
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

		label stuckToWall = 0;
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
			stuckToWall,
			wallTemperature,
			wallVectors,
			vibLevel
		);
	}
	else //Received point outside of mesh, this is almost certainly a rounding error, so snap to nearest mesh boundary to avoid loosing parcel
	{
		point newPosition = position;

		if(newPosition[0] < meshMin_[0])
		{
			newPosition[0] = meshMin_[0];
		}

		if(newPosition[0] > meshMax_[0])
		{
			newPosition[0] = meshMax_[0];
		}

		if(newPosition[1] < meshMin_[1])
		{
			newPosition[1] = meshMin_[1];
		}

		if(newPosition[1] > meshMax_[1])
		{
			newPosition[1] = meshMax_[1];
		}

		if(newPosition[2] < meshMin_[2])
		{
			newPosition[2] = meshMin_[2];
		}

		if(newPosition[2] > meshMax_[2])
		{
			newPosition[2] = meshMax_[2];
		}

		cell = -1;
		tetFace = -1;
		tetPt = -1;

		mesh_.findCellFacePt
		(
			newPosition,
			cell,
			tetFace,
			tetPt
		);

		if(cell != -1)
		{
			std::cout << "WARNING. Received parcel position snapped to mesh = ("
					  << newPosition[0] << "," << newPosition[1] << "," << newPosition[2] << ")"
					  << std::endl;

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

			label stuckToWall = 0;
			scalarField wallTemperature(4, 0.0);
			vectorField wallVectors(4, vector::zero);

			cloud_.addNewParcel
			(
				newPosition,
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
				stuckToWall,
				wallTemperature,
				wallVectors,
				vibLevel
			);
		}
		else //Position still outside of mesh (something very odd happening)
		{
			std::cout << "WARNING. Received parcel position outside of mesh and not added = ("
					  << newPosition[0] << "," << newPosition[1] << "," << newPosition[2] << ")"
					  << std::endl;
		}
	}
}

} // End namespace Foam

// ************************************************************************* //
