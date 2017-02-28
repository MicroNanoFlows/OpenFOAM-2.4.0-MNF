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
    rectangularBorder

Description
    NOT PARALLISED - massive problems here


SourceFiles
    rectangularBorderI.H
    rectangularBorder.C
    rectangularBorderIO.C

\*---------------------------------------------------------------------------*/

#ifndef rectangularBorder_H
#define rectangularBorder_H

#include "borderModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "agentCloud.H"
#include "agentWallForce.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class rectangularBorder Declaration
\*---------------------------------------------------------------------------*/

class rectangularBorder
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
        
        bool addWallForce_;
        
        autoPtr<agentWallForce> wallModel_;


        
    // Private Member Functions
        
        void initialiseBorders();

        void checkClosedEndedBorders();
       
        void reflect(label index, agent* a); 

        bool isPointWithinBorder(label index, const vector& r);       
        
public:

    //- Runtime type information
    TypeName("rectangularBorder");

    // Constructors

        //- Construct from components
        rectangularBorder
        (
            Time& time,
            agentCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~rectangularBorder();

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
