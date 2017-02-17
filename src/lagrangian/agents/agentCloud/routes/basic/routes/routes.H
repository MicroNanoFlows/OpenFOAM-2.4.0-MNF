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
    routes

Description



    

\*---------------------------------------------------------------------------*/

#ifndef routes_H
#define routes_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "routeModel.H"
#include "agent.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class routes Declaration
\*---------------------------------------------------------------------------*/

class routes
{

    // Private data


        const fvMesh& mesh_;
        
        Time& time_;
        
        agentCloud& cloud_;

//         const agentProperties& cP_;
        
        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary routesDict_;
     
        
        //- routes
        PtrList<entry> routeList_;
        List<word> routeNames_;
        List<label> routeIds_;
        List< autoPtr<routeModel> > routes_;
        
        
    // Private Member Functions

        
public:

    friend class agent;    
    
    // Constructors

        //- Constructor 
        routes
        (
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud
        );

        ~routes();

    // Member Functions

        //- initial configuration
        void initialConfig();

//         void setSchedules();
//         
//         void afterForce();
// 
//         
//         void write();
        
        // Access

//             const label& nStateControllers() const;

            inline List< autoPtr<routeModel> >& routeModels();
            inline const List< autoPtr<routeModel> >& routeModels() const;
// 
//             inline const List<word>& stateControllersNames() const; 
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "routesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
