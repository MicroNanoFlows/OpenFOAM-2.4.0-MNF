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
    schedules

Description



    

\*---------------------------------------------------------------------------*/

#ifndef schedules_H
#define schedules_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "scheduleModel.H"
#include "agent.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class schedules Declaration
\*---------------------------------------------------------------------------*/

class schedules
{

    // Private data


        const fvMesh& mesh_;
        
        Time& time_;
        
        agentCloud& cloud_;

        const agentProperties& cP_;
        
        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary schedulesDict_;
     
        
        //- schedules
        PtrList<entry> scheduleList_;
        List<word> scheduleNames_;
        List<label> scheduleIds_;
        List< autoPtr<scheduleModel> > schedules_;
        
        
    // Private Member Functions

        
public:

    friend class agent;    
    
    // Constructors

        //- Constructor 
        schedules
        (
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud
        );

        ~schedules();

    // Member Functions

        //- initial configuration
//         void initialConfig();
// 
//         void afterMove();
//         
//         void afterForce();
// 
//         
//         void write();
        
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

#include "schedulesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
