/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    newFaceZone

Description
    set a face zone 

SourceFiles
    newFaceZoneI.H
    newFaceZone.C
    newFaceZoneIO.C

\*---------------------------------------------------------------------------*/

#ifndef newFaceZone_H
#define newFaceZone_H


// #include "Time.H"
#include "IOdictionary.H"
#include "IOstreams.H"
#include "polyMesh.H"
#include "boundedBox.H"
#include "Switch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class newFaceZone Declaration
\*---------------------------------------------------------------------------*/

class newFaceZone
{
    // Private data

        const polyMesh& mesh_;

        const dictionary& dict_;
        
        word regionName_;

        labelField faces_;

        word option_;

        bool writeFaceSets_;
    
    // Private Member Functions

        
public:

    // Constructors

        //- Construct from properties

        newFaceZone
        (
            const polyMesh& mesh,
            const dictionary& dict
        );
    // Destructor

        ~newFaceZone();

    // Member Functions

        void setZone();

        // Access

        inline const labelField& faces() const;

        inline const word& name() const;

        inline const bool& writeFaceSets() const;
        // Check

        // Edit

        // Write


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "newFaceZoneI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
