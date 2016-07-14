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
    agentProperties

Description


SourceFiles
    agentProperties.C

\*---------------------------------------------------------------------------*/

#ifndef agentProperties_H
#define agentProperties_H

#include "Time.H"
#include "polyMesh.H"
#include "diagTensor.H"
#include "transform.H"
// #include "rotationTensor.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agentProperties Declaration
\*---------------------------------------------------------------------------*/

class agentProperties
{
    // Private data

        const polyMesh& mesh_;
        
        // number of agent types 
        label N_;
        
        // ordered list of agent names
        List<word> names_;
        
        // list of total masses of each agent type
        List<scalar> masses_;
        
        
        // any potential interactions go here 
        
public:

    // Constructors

        //- constructor
        agentProperties
        (
            const polyMesh& mesh
        );

    // Destructor

        ~agentProperties();


    // Member functions
        
        // no of id types
        const label& nAgentTypes() const;
        
        // namelist of all ids 
        const List<word>& agentIds() const;
        
        const List<scalar>& mass() const;
        const scalar& mass(const label& id) const;        

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
