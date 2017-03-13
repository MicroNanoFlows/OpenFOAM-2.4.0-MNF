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
    

Description

SourceFiles
    wallCompressionI.H
    wallCompression.C
    wallCompressionIO.C

\*---------------------------------------------------------------------------*/

#ifndef wallCompression_H
#define wallCompression_H


#include "agentWallForce.H"
#include "agentCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class wallCompression Declaration
\*---------------------------------------------------------------------------*/

class wallCompression
:
    public agentWallForce
{

private:

    // Private data

        dictionary propsDict_;

//         vector force_;

        scalar A_;
        scalar kappa_;
        
public:

    //- Runtime type information
    TypeName("wallCompression");

    // Static data members

    // Constructors

        //- construct from components
        wallCompression
        (
            Time& time,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~wallCompression();

    // Member Functions


        // Access

        // Check

        // Edit

            vector force(agent* agentI, const vector& rW);
            
            vector force(agent* agentI, const vector& rW, const scalar& time);
            
            vector force(agent* agentI, const scalar& time);
            
            void updateForce();
            
        // Write


            void write
            (
                const fileName& fixedPathName,
                const fileName& timePath
            );

    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
