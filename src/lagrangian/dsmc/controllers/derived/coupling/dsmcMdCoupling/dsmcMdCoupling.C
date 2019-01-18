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
    couplingInterface1d &oneDInterfaces,
    couplingInterface2d &twoDInterfaces,
    couplingInterface3d &threeDInterfaces
)
:
    dsmcCouplingController(t, cloud, dict, oneDInterfaces, twoDInterfaces, threeDInterfaces),
    propsDict_(dict.subDict(typeName + "Properties")),
    propsDictSend_(dict.subDict(typeName + "Sending")),
    propsDictRecv_(dict.subDict(typeName + "Receiving")),
    molIds_(),
    output_(false),
    oneDInterfaces_(),
    twoDInterfaces_(),
    threeDInterfaces_(threeDInterfaces),
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
    sendInterfaces_.clear();
    sendInterfaces_.setSize(0);

    if(propsDictSend_.found("sendingInterfaces"))
    {
        const List<word> interfaces(propsDictSend_.lookup("sendingInterfaces"));
        sendInterfaces_.setSize(interfaces.size(), NULL);
        sendInterfaceNames_.setSize(interfaces.size());

        for(size_t i=0; i<interfaces.size(); ++i)
        {
            //- Find MUI interfaces
            for(size_t j=0; j<threeDInterfaces.interfaces->size(); ++j)
            {
                //- If the MUI interface is found then create a copy of its pointer address and store in sendInterfaces_
                if(threeDInterfaces.interfaces->getInterfaceName(j).compare(interfaces[i]) == 0)
                {
                    sendInterfaces_[i] = threeDInterfaces.interfaces->getInterface(j);
                    sendInterfaceNames_[i] = interfaces[i]; //- Store the receiving interface name
                    break;
                }
            }
        }

        //- Check all interfaces were found
        for(size_t i=0; i<sendInterfaces_.size(); ++i)
        {
            if(sendInterfaces_[i] == NULL)
            {
                FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                            << "Could not find 3D MUI coupling interface (" << interfaces[i]
                            << ") to send for domain " << threeDInterfaces.domainName << exit(FatalError);
            }
            else
            {
                Info << "dsmcMdCoupling::dsmcMdCoupling(): Found 3D MUI coupling interface ("
                     << interfaces[i] << ") to send for domain " << threeDInterfaces.domainName << endl;
            }
        }
    }

    //- Determine receiving interfaces if defined
    recvInterfaces_.clear();
    recvInterfaces_.setSize(0);

    if(propsDictRecv_.found("receivingInterfaces"))
    {
        const List<word> interfaces(propsDictRecv_.lookup("receivingInterfaces"));
        recvInterfaces_.setSize(interfaces.size(), NULL);
        recvInterfaceNames_.setSize(interfaces.size());

        for(size_t i=0; i<interfaces.size(); ++i)
        {
            recvInterfaces_[i] = NULL;
            //- Find MUI interfaces
            for(size_t j=0; j<threeDInterfaces.interfaces->size(); ++j)
            {
                //- If the MUI interface is found then create a copy of its pointer address and store in sendInterfaces_
                if(threeDInterfaces.interfaces->getInterfaceName(j).compare(interfaces[i]) == 0)
                {
                    recvInterfaces_[i] = threeDInterfaces.interfaces->getInterface(j);
                    recvInterfaceNames_[i] = interfaces[i]; //- Store the receiving interface name
                    break;
                }
            }
        }

        //- Check all interfaces were found
        for(size_t i=0; i<recvInterfaces_.size(); ++i)
        {
            if(recvInterfaces_[i] == NULL)
            {
                FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                            << "Could not find 3D MUI coupling interface (" << interfaces[i]
                            << ") to receive for domain " << threeDInterfaces.domainName << exit(FatalError);
            }
            else
            {
                Info << "dsmcMdCoupling::dsmcMdCoupling(): Found 3D MUI coupling interface ("
                     << interfaces[i] << ") to receive for domain " << threeDInterfaces.domainName << endl;
            }
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

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    if(propsDict_.found("binModel"))
    {
        binModel_ =  autoPtr<binModel>
        (
            binModel::New(mesh_, propsDict_)
        );
    }

    if (propsDict_.found("output"))
    {
        output_ = Switch(propsDict_.lookup("output"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcMdCoupling::~dsmcMdCoupling()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcMdCoupling::initialConfiguration()
{}

void dsmcMdCoupling::controlParcelsBeforeMove()
{
    if(sending_)
    {
        sendCoupledRegion();
    }
}

void dsmcMdCoupling::calculateProperties()
{
    if(sending_)
    {
        sendCoupledParcels();
    }

    if(receiving_)
    {
        receiveCoupledMolecules();
    }
}

void dsmcMdCoupling::sendCoupledRegion()
{
#ifdef USE_MUI
    //- Only send data if at least one sending interface is defined
    if(sending_)
    {
        dsmcParcel* parcelI = NULL;

        // Iterate through all sending interfaces for this controller
        forAll(sendInterfaces_, iface)
        {
            forAll(regionIds(), id)
            {
                forAll(controlZone(regionIds()[id]), c)
                {
                    const label& cellI = controlZone(regionIds()[id])[c];
                    const List<dsmcParcel*>& parcelsInCell = cloud_.cellOccupancy()[cellI];

                    forAll(parcelsInCell, p) //Iterate through parcels in cell
                    {
                        parcelI = parcelsInCell[p];

                        mui::point3d molCentre(parcelI->position());
                        vector molVel(parcelI->U());

                        sendInterfaces_[iface]->push("mol_vel_x", molCentre, molVel[0]);
                        sendInterfaces_[iface]->push("mol_vel_y", molCentre, molVel[1]);
                        sendInterfaces_[iface]->push("mol_vel_z", molCentre, molVel[2]);
                    }
                }
            }
        }
    }
#endif
}

void dsmcMdCoupling::sendCoupledParcels()
{

}

void dsmcMdCoupling::receiveCoupledMolecules()
{

}

void dsmcMdCoupling::output
(
    const fileName& fixedPathName,
    const List<fileName>& timePaths
)
{
    /*
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {

    }
    */
}

void dsmcMdCoupling::updateProperties(const dictionary& newDict)
{}

} // End namespace Foam

// ************************************************************************* //
