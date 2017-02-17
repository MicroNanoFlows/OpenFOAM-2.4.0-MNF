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

Class
    boxInitialise

Description


SourceFiles
    boxInitialiseI.H
    boxInitialise.C
    boxInitialiseIO.C


\*---------------------------------------------------------------------------*/

#ifndef boxInitialise_H
#define boxInitialise_H

#include "agentConfiguration.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "agentCloud.H"
#include "boundedBox.H"
#include "distributePoints.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class boxInitialise Declaration
\*---------------------------------------------------------------------------*/

class boxInitialise
:
	public agentConfiguration
{

private:

    // Private data
    
        //- subDictionary containing the properties
        dictionary propsDict_;


    // Private Member Functions

        void setBoundBox
        (
            const dictionary& propsDict,
            boundedBox& bb,
            const word& name 
        );

public:

    //- Runtime type information
    TypeName("boxInitialise");

    // Constructors

        //- Construct from components
        boxInitialise
        (
            agentCloud& molCloud,
            const dictionary& dict
        );

    // Destructor

        ~boxInitialise();


    // Member Functions

        void setInitialConfiguration();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
