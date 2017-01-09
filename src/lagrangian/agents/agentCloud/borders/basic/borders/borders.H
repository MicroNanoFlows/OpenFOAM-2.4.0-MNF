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
    borders

Description



    

\*---------------------------------------------------------------------------*/

#ifndef borders_H
#define borders_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "borderModel.H"
#include "agent.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class borders Declaration
\*---------------------------------------------------------------------------*/

class borders
{

    // Private data


        const fvMesh& mesh_;
        
        Time& time_;
        
        agentCloud& cloud_;

        const agentProperties& cP_;
        
        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary bordersDict_;
     
        
        //- borders
        PtrList<entry> borderList_;
        List<word> borderNames_;
        List<label> borderIds_;
        List< autoPtr<borderModel> > borders_;
        
        
    // Private Member Functions

        
public:

    friend class agent;    
    
    // Constructors

        //- Constructor 
        borders
        (
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud
        );

        ~borders();

    // Member Functions

        //- initial configuration
        void initialConfig();

        void afterMove();
        
        void afterForce();

        
        void write();
        
        // Access

//             const label& nStateControllers() const;

//             inline List< autoPtr<pairPotential> >& stateControllers();
//             inline const List< autoPtr<pairPotential> >& stateControllers() const;
// 
//             inline const List<word>& stateControllersNames() const; 
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bordersI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
