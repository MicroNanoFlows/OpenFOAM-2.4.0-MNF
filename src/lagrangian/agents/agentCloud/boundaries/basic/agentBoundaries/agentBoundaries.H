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
    agentBoundaries

Description

	Stores all the information for the boundary models. Models are defined
    within the boundariesDict. This class selects & builds the models 
    automatically. We assume that boundaries are divided into three:

    1. Patch models -- actions on instantaneous collisions with a patch boundary 

    2. cyclic boundary models -- actions on instantaneous collisions with 
                                a cyclic boundary 

    3. General boundary models -- time-dependent actions on molecules close to
                                  the boundary (e.g. fluxes)


\*---------------------------------------------------------------------------*/

#ifndef agentBoundaries_H
#define agentBoundaries_H

#include "List.H"

#include "agentCyclicBoundary.H"
// #include "agentPatchBoundary.H"
// #include "agentGeneralBoundary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agentBoundaries Declaration
\*---------------------------------------------------------------------------*/

class agentBoundaries
{
    // Private data

        Time& time_;

        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary agentBoundariesDict_;

        label nCyclicBoundaryModels_;

        //- cyclic boundary models
        PtrList<entry> cyclicBoundaryList_;
        List<word> cyclicBoundaryNames_;
        List<label> cyclicBoundaryIds_;
        List<fileName> cMFixedPathNames_;
        List< autoPtr<agentCyclicBoundary> > cyclicBoundaryModels_;
        List<label> cyclicBoundaryToModelId_;

        
                
        void checkCyclicBoundaryModels(const polyMesh& mesh);

public:

    // Constructors
    
        //- Null Constructor 
        agentBoundaries
        (
            Time& t,
            const polyMesh& mesh
        );

        //- Constructor 
        agentBoundaries
        (
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud
        );

        ~agentBoundaries();

        void initialConfig();

//         void controlAfterForces();

        void calculateProperties();

//         void controlAfterMove();

        void outputResults();



    // Member Functions

        const label& nCyclicBoundaryModels() const;

        inline List< autoPtr<agentCyclicBoundary> >& cyclicBoundaryModels();
        inline const List< autoPtr<agentCyclicBoundary> >& cyclicBoundaryModels() const;
        inline const List<word>& cyclicBoundaryNames() const;
        inline const List<label>& cyclicBoundaryIds() const;
        inline const List<label>& cyclicBoundaryToModelIds() const;





    // Ostream Operator

//         friend Ostream& operator<<(Ostream&, const agentBoundaries&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "agentBoundariesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
