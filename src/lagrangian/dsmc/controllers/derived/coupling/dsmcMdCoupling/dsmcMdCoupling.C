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
#else
    FatalErrorIn("dsmcMdCoupling::dsmcMdCoupling()")
                << "MUI library not enabled at compilation" << exit(FatalError);
#endif

    //Determine sending properties
    if(propsDictSend_.found("mass"))
    {
        sendMass_ = Switch(propsDictSend_.lookup("mass"));
    }

    if(propsDictSend_.found("density"))
    {
        sendDensity_ = Switch(propsDictSend_.lookup("density"));
    }

    if((sendInterfaces_.size() != 0) || sendMass_ || sendDensity_)
    {
        sending_ = true;
    }

    if(propsDictRecv_.found("mass"))
    {
        recvMass_ = Switch(propsDictRecv_.lookup("mass"));
    }

    if(propsDictRecv_.found("density"))
    {
        recvDensity_ = Switch(propsDictRecv_.lookup("density"));
    }

    if((recvInterfaces_.size() != 0) || recvMass_ || recvDensity_)
    {
        receiving_ = true;

        recvMassValues_.setSize(recvInterfaces_.size());
        recvDensityValues_.setSize(recvInterfaces_.size());
    }

    writeInTimeDir_ = true;
    writeInCase_ = true;

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
{
    receiveInitialisation();
}

void dsmcMdCoupling::controlParcelsBeforeMove()
{
    sendCoupledRegion();
}

void dsmcMdCoupling::calculateProperties()
{
    sendCoupledParcels();
    receiveCoupledMolecules();
}

void dsmcMdCoupling::receiveInitialisation()
{

}

void dsmcMdCoupling::sendCoupledRegion()
{

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
    const fileName& timePath
)
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {

    }
}

void dsmcMdCoupling::updateProperties(const dictionary& newDict)
{}

} // End namespace Foam

// ************************************************************************* //
