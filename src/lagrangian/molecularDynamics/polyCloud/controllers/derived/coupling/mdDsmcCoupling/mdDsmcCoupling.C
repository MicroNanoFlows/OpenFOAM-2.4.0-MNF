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
	rU_(molCloud_.redUnits())
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
        	Info << "mdDsmcCoupling::mdDsmcCoupling(): Found 3D MUI coupling interface ("
				 << interfaces[i] << ") to send for domain " << threeDInterfaces.domainName << endl;
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
        	Info << "mdDsmcCoupling::mdDsmcCoupling(): Found 3D MUI coupling interface ("
				 << interfaces[i] << ") to receive for domain " << threeDInterfaces.domainName << endl;
        }

        molChanged_.setSize(recvInterfaces_.size());
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

    const List<word> types(propsDict_.lookup("molIds"));
	molNames_.setSize(types.size());

	forAll(types, type)
	{
		molNames_[type] = types[type];
	}

	forAll(molNames_, molType)
	{
		const label molId = findIndex(idList, molNames_[molType]);

		if(molId == -1)
		{
			FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
				<< "Cannot find molecule id: " << molNames_[molType] << nl << "in idList."
				<< exit(FatalError);
		}
	}

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

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
		//Initialise exact time sampler for MUI with a numerical tolerance of 1e-9, large value needed as dsmcFoamPlus works in non-normalised numerics.
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

mdDsmcCoupling::~mdDsmcCoupling()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mdDsmcCoupling::initialConfiguration()
{
	if(receiving_)
    {
		receiveCoupledRegion(true); // Receive ghost molecules in coupled regions at time = startTime
    }
}

void mdDsmcCoupling::controlAfterMove(int stage)
{
	if(stage == 1)
	{
		if(sending_)
		{
			sendCoupledMolecules(); // Send any molecules deleted by coupling boundary (non-blocking)
		}

		if(receiving_)
		{
			receiveCoupledParcels(); // Receive any molecules from dsmc coupling boundary (blocking)
		}
	}
	else if (stage == 2)
	{
		receiveCoupledRegion(false); // Receive ghost molecules in coupled region (blocking)
    }
}

void mdDsmcCoupling::calculateProperties()
{}

void mdDsmcCoupling::receiveCoupledRegion(bool init)
{
#ifdef USE_MUI
	label molCount = 0;
	scalar couplingTime;
	if(init)
	{
		couplingTime = 1.0;
	}
	else
	{
		couplingTime = time_.time().value() * oneOverRefTime_;
	}
	List<std::vector<mui::point3d> > rcvPoints(recvInterfaces_.size());
    List<std::vector<scalar> > rcvVelX(recvInterfaces_.size());
    List<std::vector<scalar> > rcvVelY(recvInterfaces_.size());
    List<std::vector<scalar> > rcvVelZ(recvInterfaces_.size());
    List<std::vector<scalar> > rcvMolChanged(recvInterfaces_.size());
	std::stringstream rcvStr;
	const scalar temperature = 1; //This needs to be properly defined, only correct for sample Argon case

	// Iterate through all receiving interfaces for this controller and extract a points list for each molecule type handled
	forAll(recvInterfaces_, iface)
	{
		forAll(molNames_, molType)
		{
			rcvStr.str("");
			rcvStr.clear();
			rcvStr << molNames_[molType] << "_" << "parc_changed"; //Receive string in format [type]_parc_changed

			//- Extract a list of all molecule locations received from other solver through this interface
			rcvPoints[iface] = recvInterfaces_[iface]->fetch_points<scalar>(rcvStr.str(), couplingTime, *chrono_sampler);			
                        
			//- Extract a list of all molecule change status values received from other solver through this interface
			rcvMolChanged[iface] = recvInterfaces_[iface]->fetch_values<scalar>(rcvStr.str(), couplingTime, *chrono_sampler);		

			//- Extract a list of all molecule velocities received from other solver through this interface
            rcvStr.str("");
			rcvStr.clear();
			rcvStr << molNames_[molType] << "_" << "parc_vel_x_region";
			rcvVelX[iface] = recvInterfaces_[iface]->fetch_values<scalar>(rcvStr.str(), couplingTime, *chrono_sampler);
			rcvStr.str("");
			rcvStr.clear();
			rcvStr << molNames_[molType] << "_" << "parc_vel_y_region";
			rcvVelY[iface] = recvInterfaces_[iface]->fetch_values<scalar>(rcvStr.str(), couplingTime, *chrono_sampler);
			rcvStr.str("");
			rcvStr.clear();
			rcvStr << molNames_[molType] << "_" << "parc_vel_z_region";
			rcvVelZ[iface] = recvInterfaces_[iface]->fetch_values<scalar>(rcvStr.str(), couplingTime, *chrono_sampler);			
		}
	}

	forAll(rcvPoints, ifacepts)
	{
		if(rcvPoints[ifacepts].size() > 0)
		{
			bool newList = false;		
			
			//If list size has changed then treat this as a new list
			if(static_cast<size_t>(molChanged_[ifacepts].size()) != rcvPoints[ifacepts].size())
			{
				if(molChanged_[ifacepts].size() > 0)
				{				
					forAll(molHistory_[ifacepts], mol)
					{
						if(molHistory_[ifacepts][mol] != NULL)
						{
							molCloud_.deleteParticle(*molHistory_[ifacepts][mol]);
						}
					}

					molChanged_[ifacepts].clear();
					molHistory_[ifacepts].clear();
				}
				molChanged_[ifacepts].setSize(rcvPoints[ifacepts].size(), false);
				molHistory_[ifacepts].setSize(rcvPoints[ifacepts].size(), NULL);
				newList = true;
			}

			forAll(molNames_, molType)
			{
				const label molId = findIndex(idList, molNames_[molType]);

				if(molId != -1)
				{
					vector exactBoundaryMin(vector::zero), exactBoundaryMax(vector::zero);

					if(recvInterfaceNames_[ifacepts] == "ifs_1")
					{
						exactBoundaryMin[0] = 279.411765;
						exactBoundaryMin[1] = 0;
						exactBoundaryMin[2] = 0;

						exactBoundaryMax[0] = 294.117647;
						exactBoundaryMax[1] = 147.058824;
						exactBoundaryMax[2] = 147.058824;
					}

					if(recvInterfaceNames_[ifacepts] == "ifs_2")
					{
						exactBoundaryMin[0] = 441.176471;
						exactBoundaryMin[1] = 0;
						exactBoundaryMin[2] = 0;

						exactBoundaryMax[0] = 455.882353;
						exactBoundaryMax[1] = 147.058824;
						exactBoundaryMax[2] = 147.058824;
					}

					for (size_t pts = 0; pts < rcvPoints[ifacepts].size(); pts++)
					{
						if(rcvPoints[ifacepts][pts][0] == 0.0 || rcvPoints[ifacepts][pts][1] == 0.0 || rcvPoints[ifacepts][pts][2] == 0.0)
						{
							std::cout << "receiveCoupledRegion(): [" << rcvPoints[ifacepts][pts][0] << "," << rcvPoints[ifacepts][pts][1] << "," << rcvPoints[ifacepts][pts][2] << "]" << std::endl;
						}

						molChanged_[ifacepts][pts] = static_cast<bool>(rcvMolChanged[ifacepts][pts]);
		                                
						vector velocity;
		                velocity[0] = rcvVelX[ifacepts][pts] / rU_.refVelocity();
						velocity[1] = rcvVelY[ifacepts][pts] / rU_.refVelocity();
						velocity[2] = rcvVelZ[ifacepts][pts] / rU_.refVelocity();
						
						point checkedPosition(rcvPoints[ifacepts][pts][0] * refLength_, rcvPoints[ifacepts][pts][1] * refLength_, rcvPoints[ifacepts][pts][2] * refLength_);

						if(checkedPosition[0] < exactBoundaryMin[0])
						{
							std::cout << "Trunc 0" << std::endl;
							checkedPosition[0] = exactBoundaryMin[0];
						}

						if(checkedPosition[0] > exactBoundaryMax[0])
						{
							std::cout << "Trunc 1" << std::endl;
							checkedPosition[0] = exactBoundaryMax[0];
						}

						if(checkedPosition[1] < exactBoundaryMin[1])
						{
							std::cout << "Trunc 2" << std::endl;
							checkedPosition[1] = exactBoundaryMin[1];
						}

						if(checkedPosition[1] > exactBoundaryMax[1])
						{
							std::cout << "Trunc 3" << std::endl;
							checkedPosition[1] = exactBoundaryMax[1];
						}

						if(checkedPosition[2] < exactBoundaryMin[2])
						{
							std::cout << "Trunc 4" << std::endl;
							checkedPosition[2] = exactBoundaryMin[2];
						}

						if(checkedPosition[2] > exactBoundaryMax[2])
						{
							std::cout << "Trunc 5" << std::endl;
							checkedPosition[2] = exactBoundaryMax[2];
						}

						const point position(checkedPosition[0], checkedPosition[1], checkedPosition[2]);

						if(newList) //This is a completely new list so all molecules to be inserted regardless of molChanged flag
						{
							molHistory_[ifacepts][pts] = insertMolecule(position, molId, true, temperature, velocity);
							molCount++;
						}
						else //This is not a completely new list so just check for individual molecule changes
						{
							if(molChanged_[ifacepts][pts]) //This molecule has changed in the list since the last time
							{							
								if(molHistory_[ifacepts][pts] != NULL)
								{
									molCloud_.deleteParticle(*molHistory_[ifacepts][pts]); //First delete the old molecule in this list position
								}
								molHistory_[ifacepts][pts] = insertMolecule(position, molId, true, temperature, velocity); //Insert the new molecule
								molCount++;					
							}
							else //This molecule already exists in the list so just update properties
							{
								if(molHistory_[ifacepts][pts] != NULL)
								{
									molHistory_[ifacepts][pts]->position()[0] = position[0];
									molHistory_[ifacepts][pts]->position()[1] = position[1];
									molHistory_[ifacepts][pts]->position()[2] = position[2];

									molHistory_[ifacepts][pts]->v()[0] = velocity[0];
									molHistory_[ifacepts][pts]->v()[1] = velocity[1];
									molHistory_[ifacepts][pts]->v()[2] = velocity[2];
									
									molCount++;
								}
							}
						}
					}
				}
			}			
		}
	}

	if(init)
	{
		forAll(recvInterfaces_, iface)
		{
			recvInterfaces_[iface]->commit(couplingTime);
		}
	}

	if(molCount > 0)
	{
		std::cout << "Number of molecules in coupled region = " << molCount << std::endl;
	}
#endif
}

void mdDsmcCoupling::sendCoupledMolecules()
{
#ifdef USE_MUI
	const DynamicList<polyMoleculeCloud::coupledMols>& molsToSend = molCloud_.coupledMolecules();
	scalar couplingTime = time_.time().value() * oneOverRefTime_;
	std::stringstream sendStr;

	if(molsToSend.size() != 0)
	{
		forAll(molsToSend, mols)
		{
			const label typeIndex = findIndex(molNames_, molsToSend[mols].molType);

			if(typeIndex != -1)
			{
				forAll(molsToSend[mols].sendingInterfaces, interface)
				{
					const label iface = findIndex(sendInterfaceNames_, molsToSend[mols].sendingInterfaces[interface]);

					if(iface != -1)
					{
						// Get the molecule centre
						mui::point3d molCentre;
						molCentre[0] = molsToSend[mols].mol->position()[0] * oneOverRefLength_;
						molCentre[1] = molsToSend[mols].mol->position()[1] * oneOverRefLength_;
						molCentre[2] = molsToSend[mols].mol->position()[2] * oneOverRefLength_;

						sendStr.str("");
						sendStr.clear();
						sendStr << molsToSend[mols].molType << "_" << "mol_vel_x_bound"; //Send string in format [type]_mol_vel_x_bound

						// Push the molecule velocity to the interface
						sendInterfaces_[iface]->push(sendStr.str(), molCentre, molsToSend[mols].mol->v()[0] * rU_.refVelocity());
						sendStr.str("");
						sendStr.clear();
						sendStr << molsToSend[mols].molType << "_" << "mol_vel_y_bound"; //Send string in format [type]_mol_vel_y_bound
						sendInterfaces_[iface]->push(sendStr.str(), molCentre, molsToSend[mols].mol->v()[1] * rU_.refVelocity());
						sendStr.str("");
						sendStr.clear();
						sendStr << molsToSend[mols].molType << "_" << "mol_vel_z_bound"; //Send string in format [type]_mol_vel_z_bound
						sendInterfaces_[iface]->push(sendStr.str(), molCentre, molsToSend[mols].mol->v()[2] * rU_.refVelocity());

						std::cout << "Coupling boundary molecule pushed at: [" << molCentre[0] << "," << molCentre[1] << "," << molCentre[2] << "]" << std::endl;
					}
				}
			}
		}
	}

    forAll(sendInterfaces_, iface)
    {
    	// Commit values to the coupling interface
    	sendInterfaces_[iface]->commit(couplingTime);
    }
#endif
}

void mdDsmcCoupling::receiveCoupledParcels()
{
#ifdef USE_MUI
	scalar couplingTime = time_.time().value() * oneOverRefTime_;
	List<std::vector<mui::point3d> > rcvPoints(recvInterfaces_.size());
    List<std::vector<scalar> > rcvVelX(recvInterfaces_.size());
    List<std::vector<scalar> > rcvVelY(recvInterfaces_.size());
    List<std::vector<scalar> > rcvVelZ(recvInterfaces_.size());
	const scalar temperature = 1; //This needs to be properly defined, only correct for sample Argon case
	std::stringstream rcvStr;

	// Iterate through all receiving interfaces for this controller and extract a points list for each molecule type handled
	forAll(recvInterfaces_, iface)
	{
			forAll(molNames_, molType)
			{
					rcvStr.str("");
					rcvStr.clear();
					rcvStr << molNames_[molType] << "_" << "parc_vel_x_bound";

					//- Extract a list of all molecule locations received from other solver through this interface
					rcvPoints[iface] = recvInterfaces_[iface]->fetch_points<scalar>(rcvStr.str(), couplingTime, *chrono_sampler);

					//- Extract a list of all molecule velocities received from other solver through this interface
					rcvVelX[iface] = recvInterfaces_[iface]->fetch_values<scalar>(rcvStr.str(), couplingTime, *chrono_sampler);
					rcvStr.str("");
					rcvStr.clear();
					rcvStr << molNames_[molType] << "_" << "parc_vel_y_bound";
					rcvVelY[iface] = recvInterfaces_[iface]->fetch_values<scalar>(rcvStr.str(), couplingTime, *chrono_sampler);
					rcvStr.str("");
					rcvStr.clear();
					rcvStr << molNames_[molType] << "_" << "parc_vel_z_bound";
					rcvVelZ[iface] = recvInterfaces_[iface]->fetch_values<scalar>(rcvStr.str(), couplingTime, *chrono_sampler);
			}
	}

        // Iterate through all receiving interfaces for this controller
	forAll(recvInterfaces_, ifacepts)
	{
		forAll(molNames_, molType)
		{
			if(rcvPoints.size() > 0)
			{
				const label molId = findIndex(idList, molNames_[molType]);
				
				if(molId != -1)
				{
					vector exactBoundaryMin(vector::zero), exactBoundaryMax(vector::zero);

					if(recvInterfaceNames_[ifacepts] == "ifs_1")
					{
						exactBoundaryMin[0] = 279.411765;
						exactBoundaryMin[1] = 0;
						exactBoundaryMin[2] = 0;

						exactBoundaryMax[0] = 279.411765;
						exactBoundaryMax[1] = 147.058824;
						exactBoundaryMax[2] = 147.058824;
					}

					if(recvInterfaceNames_[ifacepts] == "ifs_2")
					{
						exactBoundaryMin[0] = 441.176471;
						exactBoundaryMin[1] = 0;
						exactBoundaryMin[2] = 0;

						exactBoundaryMax[0] = 441.176471;
						exactBoundaryMax[1] = 147.058824;
						exactBoundaryMax[2] = 147.058824;
					}

					for (size_t pts = 0; pts < rcvPoints[ifacepts].size(); pts++)
					{
						vector velocity;
						velocity[0] = rcvVelX[ifacepts][pts] / rU_.refVelocity();
						velocity[1] = rcvVelY[ifacepts][pts] / rU_.refVelocity();
						velocity[2] = rcvVelZ[ifacepts][pts] / rU_.refVelocity();

						point checkedPosition(exactBoundaryMin[0], rcvPoints[ifacepts][pts][1] * refLength_, rcvPoints[ifacepts][pts][2] * refLength_);

						if(checkedPosition[1] == 0.0 || checkedPosition[2] == 0.0)
						{
							std::cout << "receiveCoupledParcels(): [" << rcvPoints[ifacepts][pts][0] * refLength_ << "," << checkedPosition[1] << "," << checkedPosition[2] << "]" << std::endl;
						}

						if(checkedPosition[1] < exactBoundaryMin[1])
						{
							checkedPosition[1] = exactBoundaryMin[1];
						}

						if(checkedPosition[1] > exactBoundaryMax[1])
						{
							checkedPosition[1] = exactBoundaryMax[1];
						}

						if(checkedPosition[2] < exactBoundaryMin[2])
						{
							checkedPosition[2] = exactBoundaryMin[2];
						}

						if(checkedPosition[2] > exactBoundaryMax[2])
						{
							checkedPosition[2] = exactBoundaryMax[2];
						}

						const point position(checkedPosition[0], checkedPosition[1], checkedPosition[2]);

						std::cout << "Coupling boundary parcel received at: [" << position[0] << "," << position[1] << "," << position[2] << "]" << std::endl;

		                insertMolecule(position, molId, false, temperature, velocity);
					}
				}
			}
		}
	}
#endif
}

void mdDsmcCoupling::output
(
    const fileName& fixedPathName,
    const List<fileName>& timePaths
)
{}

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
	scalar barrierTime = chrono_sampler->get_lower_bound(time_.time().value() * oneOverRefTime_);

	forAll(sendInterfaces_, iface)
	{
		sendInterfaces_[iface]->barrier(barrierTime);
	}
}

void mdDsmcCoupling::barrier(scalar time)
{
	forAll(sendInterfaces_, iface)
	{
		sendInterfaces_[iface]->barrier(time);
	}
}

void mdDsmcCoupling::barrier(label interface)
{
	scalar barrierTime = chrono_sampler->get_lower_bound(time_.time().value() * oneOverRefTime_);
	sendInterfaces_[interface]->barrier(barrierTime);
}

void mdDsmcCoupling::barrier(scalar time, label interface)
{
	sendInterfaces_[interface]->barrier(time);
}

void mdDsmcCoupling::forget()
{
	scalar time = time_.time().value() * oneOverRefTime_;
	forAll(recvInterfaces_, iface)
	{
		recvInterfaces_[iface]->forget(chrono_sampler->get_upper_bound(time), true);
	}
}

void mdDsmcCoupling::forget(scalar time)
{
	forAll(recvInterfaces_, iface)
	{
		recvInterfaces_[iface]->forget(chrono_sampler->get_upper_bound(time), true);
	}
}

void mdDsmcCoupling::forget(scalar time, label interface)
{
	recvInterfaces_[interface]->forget(chrono_sampler->get_upper_bound(time), true);
}

polyMolecule* mdDsmcCoupling::insertMolecule
(
    const point& position,
    const label& id,
    const bool& frozen,
    const scalar& temperature,
    vector& bulkVelocity
)
{
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
        point specialPosition(vector::zero);

        label special = 0;

        if (frozen)
        {
            specialPosition = position;

            special = polyMolecule::SPECIAL_FROZEN;
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
		    bulkVelocity,
		    vector::zero,
		    pi,
		    vector::zero,
		    specialPosition,
		    special,
		    id,
		    1.0,
		    molCloud_.getTrackingNumber()
	    );

	    molCloud_.updateNeighbouringRadii(newMol);

        return newMol;
    }
    else //Received point outside of mesh, this is almost certainly a rounding error, so snap to nearest mesh boundary to avoid loosing molecule
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
			std::cout << "WARNING. Received molecule position snapped to mesh = ("
					  << newPosition[0] << "," << newPosition[1] << "," << newPosition[2] << ")"
					  << std::endl;

			point specialPosition(vector::zero);

			label special = 0;

			if (frozen)
			{
				specialPosition = newPosition;

				special = polyMolecule::SPECIAL_FROZEN;
			}

			vector pi = vector::zero;

			tensor Q = I;

			polyMolecule* newMol = molCloud_.createMolecule
			(
				newPosition,
				cell,
				tetFace,
				tetPt,
				Q,
				bulkVelocity,
				vector::zero,
				pi,
				vector::zero,
				specialPosition,
				special,
				id,
				1.0,
				molCloud_.getTrackingNumber()
			);

			return newMol;
		}
		else //Position still outside of mesh (something very odd happening)
		{
			std::cout << "WARNING. Received molecule position outside of mesh and not added = ("
				 << newPosition[0] << "," << newPosition[1] << "," << newPosition[2] << ")"
				 << std::endl;

			 return NULL;
		}
    }
}

} // End namespace Foam

// ************************************************************************* //
