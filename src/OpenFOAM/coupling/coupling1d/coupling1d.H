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

Class
    Foam::coupling1d

Description

SourceFiles
    coupling1d.C

\*---------------------------------------------------------------------------*/

#ifndef coupling1d_H
#define coupling1d_H

#include "vector.H"
#include "scalar.H"
#include "List.H"

// NOTE:
// MUI header included if the switch -DUSE_MUI included during compilation.
#ifdef USE_MUI
    #include "mui.h"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coupling1d Declaration
\*---------------------------------------------------------------------------*/

class coupling1d
{
    // Private data
private:
    struct interfaceDetails
    {
        vector domainStart;
        vector domainEnd;
        word interfaceName;
        List<word> zoneNames;
        #ifdef USE_MUI
          mui::uniface<mui::config_1d>* mui_interface;
        #endif
        bool send;
        bool receive;
        bool smartSend;
        bool zoneExtents;
    };

    word domainName_;
    List<coupling1d::interfaceDetails> interfaces_;
    List<List<word> > zoneNames_;
    List<word> interfaceNames_;
    List<vector> domainStarts_;
    List<vector> domainEnds_;
    List<bool> send_;
    List<bool> receive_;
    List<bool> smart_send_;

public:

    // Constructors

        //- Construct without explicit domain size
        coupling1d
        (
            word domainName,
            List<List<word> >& zoneNames,
            List<word>& interfaceNames,
            List<bool>& send,
            List<bool>& receive,
            List<bool>& smart_send
        );

        //- Construct with explicit domain size
        coupling1d
        (
            word domainName,
            List<List<word> >& zoneNames,
            List<word>& interfaceNames,
            List<vector>& domainStarts,
            List<vector>& domainEnds,
            List<bool>& send,
            List<bool>& receive,
            List<bool>& smart_send
        );

    // Destructor

        ~coupling1d();

    // Member Functions
#ifdef USE_MUI
        mui::uniface<mui::config_1d>* getInterface(int index) const;
#endif
        size_t size() const;
        const vector& getInterfaceDomainStart(int index) const;
        const vector& getInterfaceDomainEnd(int index) const;
        word getInterfaceName(int index) const;
        List<word> getInterfaceZoneNames(int index) const;
        bool getInterfaceSendStatus(int index) const;
        bool getInterfaceReceiveStatus(int index) const;
        bool getInterfaceSmartSendStatus(int index) const;
        bool getInterfaceExtentsStatus(int index) const;
};

struct couplingInterface1d
{
    word domainName;
    coupling1d* interfaces;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
