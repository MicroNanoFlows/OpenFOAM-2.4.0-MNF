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

// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * //

Foam::coupling2d::coupling2d(word domainName, List<word>& interfaceNames, List<vector>& domainStarts,
                             List<vector>& domainEnds, List<bool>& send, List<bool>& receive)
{
    //
    //Ensure parameters for each interface passed to constructor
    if(interfaceNames.size() == domainStarts.size() == domainEnds.size() == send.size() == receive.size())
    {
        domainName_ = domainName;

        interfaceDetails newInterface;
        interfaces.setSize(interfaceNames.size());

        for(size_t i=0; i<interfaceNames.size(); i++)
        {
            std::vector<std::string> interfaceList;

            newInterface.interfaceName = interfaceNames[i];
            interfaceList.push_back(newInterface.interfaceName); //Need std::vector copy for MUI create_uniface function
            newInterface.domainStart = domainStarts[i];
            newInterface.domainEnd = domainEnds[i];
            newInterface.send = send[i];
            newInterface.receive = receive[i];

            #ifdef USE_MUI
              std::vector<mui::uniface<mui::config_2d>*> returnInterfaces;
              returnInterfaces = mui::create_uniface<mui::config_2d>(static_cast<std::string>(domainName_), interfaceList);
              newInterface.mui_interface = returnInterfaces[0];
            #endif

            interfaces[i] = newInterface;
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupling2d::~coupling2d()
{
    #ifdef USE_MUI
        for(size_t i=0; i<interfaces.size(); ++i)
        {
            delete interfaces[i].mui_interface;
        }
    #endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#ifdef USE_MUI
mui::uniface<mui::config_2d>* Foam::coupling2d::getInterface(int index) const
{
    return interfaces[index].mui_interface;
}
#endif

size_t Foam::coupling2d::size() const
{
    return interfaces.size();
}

const Foam::word Foam::coupling2d::getInterfaceName(int index) const
{
    return interfaces[index].interfaceName;
}

const Foam::vector& Foam::coupling2d::getInterfaceDomainStart(int index) const
{
    return interfaces[index].domainStart;
}

const Foam::vector& Foam::coupling2d::getInterfaceDomainEnd(int index) const
{
    return interfaces[index].domainEnd;
}

const bool Foam::coupling2d::getInterfaceSendStatus(int index) const
{
    return interfaces[index].send;
}

const bool Foam::coupling2d::getInterfaceReceiveStatus(int index) const
{
    return interfaces[index].receive;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
