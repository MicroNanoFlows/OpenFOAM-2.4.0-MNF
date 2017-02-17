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
    agentCount


Description


SourceFiles
    agentCountI.H
    agentCount.C
    agentCountIO.C

\*---------------------------------------------------------------------------*/

#ifndef agentCount_H
#define agentCount_H

#include "agentMeasurement.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "agentCloud.H"
#include "distribution.H"
#include "tiltedCylinder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agentCount Declaration
\*---------------------------------------------------------------------------*/

class agentCount
:
	public agentMeasurement
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        List<label> agentIds_;

        word fieldName_;

        tiltedCylinder zone_;

        DynamicList<label> tNs_;
        DynamicList<scalar> times_;
        
//         scalar nAgentsInstant_;
        scalar nAgentsCumul_;
        
        DynamicList<scalar> nAgentsCumulField_;
        DynamicList<scalar> timeField_;
        
        DynamicList<scalar> nAgentsInstantField_;
        
        
        
    // Private Member Functions

            
public:

    //- Runtime type information
    TypeName("agentCount");

    // Constructors

        //- Construct from components
        agentCount
		(
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud,
            const dictionary& dict
		);

    // Destructor

 		~agentCount();

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
