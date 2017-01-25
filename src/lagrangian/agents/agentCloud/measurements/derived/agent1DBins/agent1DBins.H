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
    agent1DBins


Description


SourceFiles
    agent1DBinsI.H
    agent1DBins.C
    agent1DBinsIO.C

\*---------------------------------------------------------------------------*/

#ifndef agent1DBins_H
#define agent1DBins_H

#include "agentMeasurement.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "agentCloud.H"
#include "distribution.H"
#include "boundedBox.H"
#include "binModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agent1DBins Declaration
\*---------------------------------------------------------------------------*/

class agent1DBins
:
	public agentMeasurement
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        List<label> agentIds_;

        word fieldName_;

        autoPtr<binModel> binModel_;        
        
        label nBins_;

//         scalar totalArea_;
        
        scalarField N_;
        scalarField mass_;
        vectorField mom_;
        vectorField vel_;
        scalarField kE_;  
        
        scalarField NField_;
        scalarField rhoNField_;
        scalarField rhoMField_;
        vectorField UField_;
        vectorField momField_;
        scalarField kEField_;
        
        scalar nTimeSteps_;

        bool resetAtOutput_;  
        
        scalar dZ_;
        
    // Private Member Functions

            
public:

    //- Runtime type information
    TypeName("agent1DBins");

    // Constructors

        //- Construct from components
        agent1DBins
		(
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud,
            const dictionary& dict
		);

    // Destructor

 		~agent1DBins();

    // Member Functions

        //- create initial field
        void createField();

        //- calculate field
        void calculateField();

        void measureDuringForceComputation
        (
            agent* molI,
            agent* molJ
        );

        void measureDuringForceComputationSite
        (
            agent* molI,
            agent* molJ,
            label sI,
            label sJ
        );

//         const propertyField& fields() const;

        //- write field
        void writeField();
}; 		


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
