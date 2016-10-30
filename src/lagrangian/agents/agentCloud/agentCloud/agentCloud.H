/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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
    Foam::agentCloud

Description

SourceFiles
    agentCloudI.H
    agentCloud.C


\*---------------------------------------------------------------------------*/

#ifndef agentCloud_H
#define agentCloud_H

#include "Cloud.H"
#include "agent.H"
#include "IOdictionary.H"
#include "trackingNumber.H"
// #include "cachedRandomMD.H"
#include "agentIntegrators.H"
#include "cyclicBoundaries.H"
#include "agentProperties.H"
#include "reducedUnits.H"
#include "agentMeasurements.H"
#include "agentControllers.H"
#include "boundedBox.H"
#include "forceFields.H"
#include "borders.H"
#include "agentBoundaries.H"
#include "schedules.H"
#include "Random.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class agentCloud Declaration
\*---------------------------------------------------------------------------*/

class agentCloud
:
    public Cloud<agent>
{

private:

    // Private data

        const polyMesh& mesh_;
        
        const agentProperties& cP_;
        
//         cachedRandomMD& rndGen_;
        Random rndGen_;
        
        agentIntegrators int_;
        
        reducedUnits rU_;
        
        List<DynamicList<agent*> > cellOccupancy_;
       
        agentMeasurements fields_; 
        
        agentControllers controllers_;
        
        trackingNumber agentTracking_;
        
        cyclicBoundaries cyclics_;
        
        forceFields f_;
                
        borders b_;

        agentBoundaries ob_;
        
        schedules s_;
        


    // Private Member Functions

        void checkMoleculesInMesh();

        //- Disallow default bitwise copy construct
        agentCloud(const agentCloud&);
        
        //- Disallow default bitwise assignment
        void operator=(const agentCloud&);

public:

    // Constructors

        //- Construct for agentFoam
        agentCloud
        (
            Time& t,
            const polyMesh& mesh,
            const word& cloudName,         
            const agentProperties& cP
        );

        //- General constructor
        agentCloud
        (
            Time& t,
            const polyMesh& mesh,
            const word& cloudName,
            const agentProperties& cP,
            const word& option,
            const bool& clearFields
        );
        
    // Static constructors

        static autoPtr<agentCloud> New
        (
            Time& t,
            const polyMesh& mesh,
            const word& cloudName,
            const agentProperties& cP
        );

        static autoPtr<agentCloud> New
        (
            Time& t,
            const polyMesh& mesh,
            const word& cloudName,
            const agentProperties& cP,
            const word& option,
            const bool& clearFields
        );


    // Member Functions

        void createAgent
        (
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const vector& v,
            const vector& d,
            const vector& f,
            const vector& specialPosition,
            const scalar& mass,
            const scalar& potentialEnergy,
            const scalar& R,
            const scalar& frac,   
            const label special,
            const label id,
            const label trackingNumber
        );
        
      
        
        //- Evolve the molecules (move, calculate forces, etc)
        
        void evolve();
        
        void move();
        

        // Access
        
            inline const polyMesh& mesh() const;        
            inline const List<DynamicList<agent*> >& cellOccupancy() const;            
            inline const agentProperties& cP() const;            
            inline Random& rndGen();            
            inline trackingNumber& agentTracking();            
//             inline clockTimer& clock();            
            
            inline const agentMeasurements& fields() const;  
            inline agentMeasurements& fields();            
            inline const reducedUnits& redUnits() const;
            
            inline const agentControllers& controllers() const; 
            inline agentControllers& controllers();
            
            inline const cyclicBoundaries& cyclics() const;

            inline forceFields& f();
            inline borders& b();            
            
            inline agentBoundaries& ob();             
            
            inline const schedules& s() const;     
            inline schedules& s();     
            
    // Member Operators

        label nAgents() const;
    
        void updateTrackingNumbersAfterRead();

        label getTrackingNumber();
        
        void resetTrackingNumbers();    
    
        void buildCellOccupancy();

        void insertMolInCellOccupancy
        (
            agent* mol
        );

        void removeMolFromCellOccupancy
        (
            const label& cellMolId,
            const label& cell
        );

        void removeMolFromCellOccupancy
        (
            agent* molI
        );

/*        void updateNeighbouringRadii
        (
            agent* molI
        ); */   
    

        
        void readNewField();    

        //- Write agent sites in XYZ format
//         void writeXYZ(const fileName& fName) const;
        
        void writeXMOL();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "agentCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
