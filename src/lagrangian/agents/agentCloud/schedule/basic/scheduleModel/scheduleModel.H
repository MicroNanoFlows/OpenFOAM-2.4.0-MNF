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
    scheduleModel

Description
    Basic/abstract class of a state controller


SourceFiles
    scheduleModelI.H
    scheduleModel.C
    scheduleModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef scheduleModel_H
#define scheduleModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "agent.H"
#include "writeTimeData.H"
#include "selectAgentIds.H"
#include "reducedUnits.H"
#include "vectorList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class scheduleModel Declaration
\*---------------------------------------------------------------------------*/

class scheduleModel
{

protected:

    // Protected data

        const fvMesh& mesh_;

        Time& time_;        
        
        agentCloud& cloud_;

        const reducedUnits& rU_;        
        
        
    // Private Member Functions


        
public:

    //- Runtime type information
    TypeName("scheduleModel");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            scheduleModel,
            dictionary,
            (
                Time& time,
                agentCloud& cloud,
                const dictionary& dict
            ),
            (time, cloud, dict)
        );
   
    // Constructors

        //- Construct from components
        scheduleModel
        (
            Time& time,
            agentCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<scheduleModel> New
        (
            Time& time,
            agentCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~scheduleModel();


    // Member Functions

        
        
        
       virtual void initialConfiguration() = 0;
        
       virtual void setSchedule() = 0;
        
//         virtual void afterForce() = 0;
        
//         const List<vectorList>& borderList() const;
//         
//         bool isPointWithinBorder(label index, const vector& r);
        
//         void writeBorders(const fileName& pathName);
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
