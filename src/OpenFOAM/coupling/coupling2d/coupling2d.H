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
    Foam::coupling2d

Description

SourceFiles
    coupling2d.C

\*---------------------------------------------------------------------------*/

#ifndef coupling2d_H
#define coupling2d_H

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
                         Class coupling2d Declaration
\*---------------------------------------------------------------------------*/

class coupling2d
{
    // Private data
private:
    struct interfaceDetails
    {
        word interfaceName;
        vector domainStart;
        vector domainEnd;
        bool send;
        bool receive;
        #ifdef USE_MUI
          mui::uniface<mui::config_2d>* mui_interface;
        #endif
    };

    word domainName_;
    List<coupling2d::interfaceDetails> interfaces;

public:

    // Constructors

        //- Construct
        coupling2d
        (
            word domainName,
            List<word>& interfaceNames,
            List<vector>& domainStarts,
            List<vector>& domainEnds,
            List<bool>& send,
            List<bool>& receive
        );

    // Destructor

        ~coupling2d();

    // Member Functions
#ifdef USE_MUI
        mui::uniface<mui::config_2d>* getInterface(int index) const;
#endif
        size_t size() const;
        const word getInterfaceName(int index) const;
        const vector& getInterfaceDomainStart(int index) const;
        const vector& getInterfaceDomainEnd(int index) const;
        const bool getInterfaceSendStatus(int index) const;
        const bool getInterfaceReceiveStatus(int index) const;
};

struct couplingInterface2d
{
    word domainName;
    coupling2d* interfaces;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
