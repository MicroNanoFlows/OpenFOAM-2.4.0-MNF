/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "coupling3d.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//Construct without explicit domain sizes
Foam::coupling3d::coupling3d
(
    word domainName,
    List<List<word> >& zoneNames,
    List<word>& interfaceNames,
    List<bool>& send,
    List<bool>& receive,
    List<bool>& smart_send
)
:
    domainName_(domainName),
    zoneNames_(zoneNames),
    interfaceNames_(interfaceNames),
    send_(send),
    receive_(receive),
    smart_send_(smart_send)
{
    interfaceDetails newInterface;
    interfaces_.setSize(interfaceNames_.size());

    forAll(interfaceNames_, i)
    {
        std::vector<std::string> interfaceList;

        newInterface.interfaceName = interfaceNames_[i];
        interfaceList.push_back(newInterface.interfaceName); //Need std::vector copy for MUI create_uniface function
        newInterface.zoneNames = zoneNames_[i];
        newInterface.send = send_[i];
        newInterface.receive = receive_[i];
        newInterface.smartSend = smart_send_[i];
        newInterface.zoneExtents = true;

        #ifdef USE_MUI
          auto returnInterfaces_ = mui::create_uniface<mui::config_3d>(static_cast<std::string>(domainName_), interfaceList);
          newInterface.mui_interface = returnInterfaces_[0].release();
        #endif

        interfaces_[i] = newInterface;
    }
}

//Construct with explicit domain sizes
Foam::coupling3d::coupling3d
(
    word domainName,
    List<List<word> >& zoneNames,
    List<word>& interfaceNames,
    List<vector>& domainStarts,
    List<vector>& domainEnds,
    List<bool>& send,
    List<bool>& receive,
    List<bool>& smart_send
)
:
    domainName_(domainName),
    zoneNames_(zoneNames),
    interfaceNames_(interfaceNames),
    domainStarts_(domainStarts),
    domainEnds_(domainEnds),
    send_(send),
    receive_(receive),
    smart_send_(smart_send)
{
    interfaceDetails newInterface;
    interfaces_.setSize(interfaceNames.size());

    forAll(interfaceNames_, i)
    {
        std::vector<std::string> interfaceList;

        newInterface.interfaceName = interfaceNames_[i];
        interfaceList.push_back(newInterface.interfaceName); //Need std::vector copy for MUI create_uniface function
        newInterface.zoneNames = zoneNames_[i];
        newInterface.domainStart = domainStarts_[i];
        newInterface.domainEnd = domainEnds_[i];
        newInterface.send = send_[i];
        newInterface.receive = receive_[i];
        newInterface.smartSend = smart_send_[i];
        newInterface.zoneExtents = false;

        #ifdef USE_MUI
          auto returnInterfaces_ = mui::create_uniface<mui::config_3d>(static_cast<std::string>(domainName_), interfaceList);
          newInterface.mui_interface = returnInterfaces_[0].release();
        #endif

        interfaces_[i] = newInterface;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupling3d::~coupling3d()
{
    #ifdef USE_MUI
        forAll(interfaces_, iface)
        {
            delete interfaces_[iface].mui_interface;
        }
    #endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#ifdef USE_MUI
mui::uniface<mui::config_3d>* Foam::coupling3d::getInterface(int index) const
{
    return interfaces_[index].mui_interface;
}
#endif

size_t Foam::coupling3d::size() const
{
    return interfaces_.size();
}

const Foam::vector& Foam::coupling3d::getInterfaceDomainStart(int index) const
{
    return interfaces_[index].domainStart;
}

const Foam::vector& Foam::coupling3d::getInterfaceDomainEnd(int index) const
{
    return interfaces_[index].domainEnd;
}

Foam::word Foam::coupling3d::getInterfaceName(int index) const
{
    return interfaces_[index].interfaceName;
}

Foam::List<Foam::word> Foam::coupling3d::getInterfaceZoneNames(int index) const
{
    return interfaces_[index].zoneNames;
}

bool Foam::coupling3d::getInterfaceSendStatus(int index) const
{
    return interfaces_[index].send;
}

bool Foam::coupling3d::getInterfaceReceiveStatus(int index) const
{
    return interfaces_[index].receive;
}

bool Foam::coupling3d::getInterfaceSmartSendStatus(int index) const
{
    return interfaces_[index].smartSend;
}

bool Foam::coupling3d::getInterfaceExtentsStatus(int index) const
{
    return interfaces_[index].zoneExtents;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
