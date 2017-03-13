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
    forceFields

Description



    

\*---------------------------------------------------------------------------*/

#ifndef forceFields_H
#define forceFields_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "agent.H"
#include "fvMesh.H"
#include "twoDimCellInteractions.H"

#include "bodyForce.H"
#include "pairPotential.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class forceFields Declaration
\*---------------------------------------------------------------------------*/

class forceFields
{

    // Private data


        const fvMesh& mesh_;
        
        Time& time_;
        
        agentCloud& cloud_;

        const agentProperties& cP_;
        
        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary forceFieldsDict_;

        fileName pathName_;        
        
        scalar rCut_;

        // number of one-body forceFields
        label nBodyForces_;
        
        //- one-body force field
        PtrList<entry> bodyForceList_;
        List<word> bodyForceNames_;
        List<label> bodyForceIds_;
        List< autoPtr<bodyForce> > bodyForces_;
        

        // number of pair forceFields
        label nPairPotentials_;
        
        //- pairPotentials
        PtrList<entry> pairPotList_;
        List<word> pairPotNames_;
        List<label> pairPotIds_;
        List< autoPtr<pairPotential> > pairPotentials_;
        
        bool onePotential_;
        
        // pair potential option 
        // (e.g. onePotential, table etc)
        word option_;
        
        
        // use this when given two molecules to compute forces/energies
        // e.g. index = pairPotentialLookUp_[idA][idB]
        //      pairPotentials_[index]
        // outer list = size of idlist
        // inner list = size of idlist 
        List<List<label> > pairPotentialLookUp_;
        
        // cellListAlgorithm
        twoDimCellInteractions<agent> iL_;
        
        //- interaction particle list
//         List<DynamicList<agent*> > ipl_;        
        
        List<vector> agentSocialForces_;

        List<vector> agentBodyForces_;
        
        List<vector> agentWallForces_;
        
    // Private Member Functions

        void testPairPotentials();
        
        inline void evaluatePair
        (
            agent* molI,
            agent* molJ
        );
        
public:

    friend class agent;    
    
    // Constructors

        //- Constructor 
        forceFields
        (
            Time& t,
            const polyMesh& mesh,
            agentCloud& cloud
        );

        ~forceFields();

    // Member Functions

        //- initial configuration
        void initialConfig();

        void calculateBodyForces();
        
        void calculatePairForces();
        
        void calculateWallForces();
        
        // Access

        inline List< autoPtr<bodyForce> >& bodyForces();

        inline const scalar& rCut() const;
        
        inline const List<vector>& agentSocialForces() const;

        inline const List<vector>& agentBodyForces() const;
        
        inline const List<vector>& agentWallForces() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "forceFieldsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
