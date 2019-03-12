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

#include "coupling2d.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupling2d::coupling2d
(
    word domainName,
    List<word>& interfaceNames,
    List<bool>& send,
    List<bool>& receive,
    List<bool>& smart_send,
    List<vector>& dom_send_start,
    List<vector>& dom_send_end,
    List<vector>& dom_rcv_start,
    List<vector>& dom_rcv_end
)
:
    domainName_(domainName),
    interfaceNames_(interfaceNames),
    send_(send),
    receive_(receive),
    smart_send_(smart_send),
    dom_send_start_(dom_send_start),
    dom_send_end_(dom_send_end),
    dom_rcv_start_(dom_rcv_start),
    dom_rcv_end_(dom_rcv_end)
{
    interfaces_.setSize(interfaceNames.size());

    forAll(interfaceNames_, i)
    {
        interfaceDetails newInterface;
        std::vector<std::string> interfaceList;

        newInterface.interfaceName = interfaceNames_[i];
        interfaceList.emplace_back(newInterface.interfaceName); //Need std::vector copy for MUI create_uniface function
        newInterface.send = send_[i];
        newInterface.receive = receive_[i];
        newInterface.smartSend = smart_send_[i];
        newInterface.dom_send_start = dom_send_start_[i];
        newInterface.dom_send_end = dom_send_end_[i];
        newInterface.dom_rcv_start = dom_rcv_start_[i];
        newInterface.dom_rcv_end = dom_rcv_end_[i];

        #ifdef USE_MUI
            auto returnInterfaces = mui::create_uniface<mui::config_2d>(static_cast<std::string>(domainName_), interfaceList);
            newInterface.mui_interface = returnInterfaces[0].release();
        #endif

        interfaces_[i] = newInterface;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupling2d::~coupling2d()
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
mui::uniface<mui::config_2d>* Foam::coupling2d::getInterface(int index) const
{
    return interfaces_[index].mui_interface;
}
#endif

size_t Foam::coupling2d::size() const
{
    return interfaces_.size();
}

Foam::word Foam::coupling2d::getInterfaceName(int index) const
{
    return interfaces_[index].interfaceName;
}

bool Foam::coupling2d::getInterfaceSendStatus(int index) const
{
    return interfaces_[index].send;
}

bool Foam::coupling2d::getInterfaceReceiveStatus(int index) const
{
    return interfaces_[index].receive;
}

bool Foam::coupling2d::getInterfaceSmartSendStatus(int index) const
{
    return interfaces_[index].smartSend;
}

Foam::vector Foam::coupling2d::getInterfaceSendDomStart(int index) const
{
    return interfaces_[index].dom_send_start;
}

Foam::vector Foam::coupling2d::getInterfaceSendDomEnd(int index) const
{
    return interfaces_[index].dom_send_end;
}

Foam::vector Foam::coupling2d::getInterfaceReceiveDomStart(int index) const
{
    return interfaces_[index].dom_rcv_start;
}

Foam::vector Foam::coupling2d::getInterfaceReceiveDomEnd(int index) const
{
    return interfaces_[index].dom_rcv_end;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
