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
    agentControllers

Description

	Stores all the information for the agentControllers models defined within 
    the agentControllersDict, and selects & builds the models automatically.

    

\*---------------------------------------------------------------------------*/

#ifndef agentControllers_H
#define agentControllers_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "agentController.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agentControllers Declaration
\*---------------------------------------------------------------------------*/

class agentControllers
{

    // Private data

        Time& time_;

        const polyMesh& mesh_;

        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary agentControllersDict_;

        // number of agentControllers
        label nStateControllers_;

        //- agentControllers
        PtrList<entry> stateControllersList_;
        List<word> sCNames_;
        List<label> sCIds_;
        List<fileName> sCFixedPathNames_;
        List< autoPtr<agentController> > stateControllers_;


        labelList controllersDuringForceComp_;
        
public:

    // Constructors

        //- Null Constructor 
        agentControllers
        (
            Time& t,
            const polyMesh& mesh
        );

        //- Constructor 
        agentControllers
        (
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud
        );


        ~agentControllers();

    // Member Functions

        //- initial configuration
        void initialConfig();

        void controlVelocitiesI();
        void controlBeforeMove();
        void controlBeforeForces();
        void controlDuringForceComputation
        (
            agent* molI,
            agent* molJ
        );
        void controlAfterForces();
        void controlVelocitiesII();
        
        
        //- calculate controller properties
        void calculateStateProps();

        //- output controller properties
        void outputStateResults();
        
        // Access

            const label& nStateControllers() const;

            inline List< autoPtr<agentController> >& stateControllers();
            inline const List< autoPtr<agentController> >& stateControllers() const;

            inline const List<word>& stateControllersNames() const; 
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "agentControllersI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
