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
    reflectiveRectangularBorderModified

Description
    NOT PARALLISED - massive problems here


SourceFiles
    reflectiveRectangularBorderModifiedI.H
    reflectiveRectangularBorderModified.C
    reflectiveRectangularBorderModifiedIO.C

\*---------------------------------------------------------------------------*/

#ifndef reflectiveRectangularBorderModified_H
#define reflectiveRectangularBorderModified_H

#include "borderModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "agentCloud.H"
// #include "agentWallForce.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class reflectiveRectangularBorderModified Declaration
\*---------------------------------------------------------------------------*/

class reflectiveRectangularBorderModified
:
	public borderModel
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        List<vectorList> borderList_;
        
        // border properties
        
        // normal vectors point outwards 
        List<vectorList> normalVectors_;
        List<vectorList> midPoints_;
        
        scalar treshold_;
 
        // outer list = nCells on mesh
        // inner list = index to borderList_ to interact with
        List<DynamicList<label> > interactionList_;  

        label scalingValue_;

        List<label> agentIds_;        
        
        // adding a wall force

//         autoPtr<agentWallForce> wallModel_;
        
//         List<vectorList> pointsOnWall_;
//         
//         bool addWallForce_;
//         
//         scalar dX_;


        
    // Private Member Functions
        
        void initialiseBorders();

        void checkClosedEndedBorders();
       
        void reflect(label index, agent* a); 

        bool isPointWithinBorder(label index, const vector& r, const scalar& agentRadius);       
        
public:

    //- Runtime type information
    TypeName("reflectiveRectangularBorderModified");

    // Constructors

        //- Construct from components
        reflectiveRectangularBorderModified
        (
            Time& time,
            agentCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~reflectiveRectangularBorderModified();

    // Member Functions
        
        void initialConfiguration();
        
        void afterMove();
        
        void afterForce();
        
                
        // write
        
        void write(const fileName& pathName);
 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
