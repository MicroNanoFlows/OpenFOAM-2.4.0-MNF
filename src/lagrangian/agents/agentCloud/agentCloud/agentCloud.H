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
#include "cachedRandomMD.H"
#include "clockTimer.H"

//#include "labelVector.H"
//#include "fileName.H"
//#include "reducedUnits.H"
//#include "polyFieldProperties.H"
//#include "polyBoundaries.H"
//#include "polyControllers.H"
//#include "polyFaceTracker.H"
//#include "cyclicBoundaries.H"
//#include "cellInteractions.H"
//#include "potentials.H" 
//#include "agentProperties.H"
//#include "polyIntegrators.H"


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
        
        cachedRandomMD& rndGen_;

//         polyIntegrators int_;
        
//         potentials p_;
        
        List<DynamicList<agent*> > cellOccupancy_;
       
        //- measurement fields
//         polyFieldProperties fields_; 
        
        //- controllers
//         polyControllers controllers_;
        
        
        trackingNumber agentTracking_;
        
       
//         cellInteractions<agent> iL_;
        
        //- interaction particle list
//         List<DynamicList<agent*> > ipl_;
        
        clockTimer clock_;


    // Private Member Functions

        void calculatePairForces();
//         void checkForOverlaps();
        void checkMoleculesInMesh();
        void setMoleculesForDirectInteraction();
        void computePairForces();
        void accelerationUpdate();
        void velocityUpdate(const scalar& trackTime);


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
            const agentProperties& cP,
            cachedRandomMD& rndGen
        );

        //- General constructor
        agentCloud
        (
            Time& t,
            const polyMesh& mesh,
            const agentProperties& cP,
            cachedRandomMD& rndGen,
            const word& option,
            const bool& clearFields
        );
        
    // Static constructors

        static autoPtr<agentCloud> New
        (
            Time& t,
            const polyMesh& mesh,
            const agentProperties& cP,
            cachedRandomMD& rndGen
        );

        static autoPtr<agentCloud> New
        (
            Time& t,
            const polyMesh& mesh,
            const agentProperties& cP,
            cachedRandomMD& rndGen,
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
            const vector& a,
            const vector& f,
            const vector& specialPosition,
            const scalar& R,
            const scalar& frac,   
            const label special,
            const label id,
            const label trackingNumber
        );
        
      
        
        //- Evolve the molecules (move, calculate forces, etc)
        
        void evolve();
        /*
        void evolveBeforeForces();
        void evolveAfterForces();
        
        void controlBeforeVelocity();
        void updateVelocity();
        void controlBeforeMove();
        void move();
        void move(const scalar& trackTime);        
        void updateAfterMove(const scalar& trackTime);        
        void controlAfterMove();
        void controlBeforeForces();
        void clearLagrangianFields();
        void calculateForce();
        void updateAcceleration();
        void controlAfterForces();
        void controlAfterVelocity();
        void postTimeStep();*/

        /*
        void prepareInteractions();*/
        /*
        void buildCellOccupancy();
        void rebuildCellOccupancy();

        inline void evaluatePair
        (
            agent* molI,
            agent* molJ
        );

        inline bool evaluatePotentialLimit
        (
            agent* molI,
            agent* molJ,
            const scalar& potentialEnergyLimit
        ) ;


        void updateRadii();
        void readNewField();
        void setIPL();
        void writeReferredCloud();
        */
        // Access
        
            inline const polyMesh& mesh() const;        
            inline const List<DynamicList<agent*> >& cellOccupancy() const;            
            inline const agentProperties& cP() const;            
            inline cachedRandomMD& rndGen();            
            inline trackingNumber& agentTracking();            
            inline clockTimer& clock();            
/*
            label nSites() const;
            inline const potentials& pot() const;
            inline potentials& pot();            
//             inline const potential& pot() const;

//             inline const List<agent::constantProperties> constProps() const;
//             inline const agent::constantProperties& constProps(label id) const;

            inline const reducedUnits& redUnits() const;
            inline const polyFieldProperties& fields() const;  
            inline polyFieldProperties& fields();
            inline const cyclicBoundaries& cyclics() const;
            inline const polyBoundaries& boundaries() const;
            inline polyBoundaries& boundaries();
            inline const polyControllers& controllersInfo() const; 
            inline polyControllers& controllersInfo();
            inline const polyFaceTracker& tracker() const;
            inline polyFaceTracker& tracker();

            inline cellInteractions<agent>& il();

            */
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
        void writeXYZ(const fileName& fName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "agentCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
