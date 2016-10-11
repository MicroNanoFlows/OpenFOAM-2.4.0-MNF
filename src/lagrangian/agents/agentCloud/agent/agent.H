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
    Foam::agent

Description
    Foam::agent

SourceFiles
    agentI.H
    agent.C
    agentIO.C

\*---------------------------------------------------------------------------*/

#ifndef agent_H
#define agent_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "agentProperties.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class agentCloud;

/*---------------------------------------------------------------------------*\
                           Class agent Declaration
\*---------------------------------------------------------------------------*/

class agent
:
    public particle
{

public:

    // Values of special that are less than zero are for built-in functionality.
    // Values greater than zero are user specifiable/expandable (i.e. test
    // special_ >= SPECIAL_USER)

    enum specialTypes
    {
        SPECIAL_TETHERED = -1,
        SPECIAL_FROZEN   = -2,
        NOT_SPECIAL      = 0,
        SPECIAL_REFERRED = 1,
        SPECIAL_USER     = 2
    };

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::TrackingData<agentCloud>
    {
        // label specifying which part of the integration algorithm is taking
        label part_;


    public:

        // Constructors

            trackingData(agentCloud& cloud, label part)
            :
                particle::TrackingData<agentCloud>(cloud),
                part_(part)
            {}

        // Member functions

            inline label part() const
            {
                return part_;
            }
    };    
    

private:

    // Private data

        //- Be careful with the ordering of data.
        //  It has an impact on binary transfer:
        //    -# Put the largest data members 1st
        //    -# Pair up labels,
        //    -# Don't go scalar-label, scalar-label, because in 64bit mode,
        //       the labels will be padded by 4bytes.

        //- linear velocity
        vector v_;

        //- destination
        vector d_;
        
        //- linear acceleration
//         vector a_;
        
        //- force
        vector f_;

       
        vector specialPosition_;

        scalar mass_;
        
        scalar potentialEnergy_;

        scalar R_;

        // - fraction of force potential contribution
        scalar frac_;    
        
        label special_;

        label id_;

        //- each agent will have a unique label
        label trackingNumber_;


    // Private Member Functions



public:

    friend class Cloud<agent>;

    // Constructors

        //- Construct from components
        inline agent
        (
            const polyMesh& mesh,            
            const agentProperties& cP,
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

        //- Construct from Istream
        agent
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new agent(*this));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<agent> operator()(Istream& is) const
            {
                return autoPtr<agent>(new agent(mesh_, is, true));
            }
        };

    // Member Functions

        // Tracking

            bool move(trackingData&, const scalar& trackTime);

            void transformProperties(const tensor& T);

            void transformProperties(const vector& separation);

//             void setSitePositions(const agentProperties& cP);

//             void setSiteSizes(label size);

//             void updateHalfVelocity
//             (
//                 const agentProperties& cP,
//                 const scalar& trackTime
//             );

//             void updateAcceleration
//             (
//                 const agentProperties& cP
//             );

//             void updateAfterMove
//             (
//                 const agentProperties& cP,
//                 const scalar& trackTime
//             );

            void setAsReferred();
            
        // Access

            inline const vector& v() const;
            inline vector& v();

            inline const vector& d() const;
            inline vector& d();            
            
            inline const vector& f() const;
            inline vector& f();            
            
            inline vector a();
            
            inline const vector& specialPosition() const;
            inline vector& specialPosition();

            inline const scalar& mass() const;
            inline scalar& mass();             
            
            inline const scalar& potentialEnergy() const;
            inline scalar& potentialEnergy();            
            
            // - return shortest radius
            inline const scalar& R() const;
            inline scalar& R();

            inline const scalar& fraction() const;
            inline scalar& fraction();

            inline label special() const;
            inline label& special();
            
            inline label id() const;
            inline label& id();            
            
            // - return tracking number
            inline const label& trackingNumber() const;
            inline label& trackingNumber();

            inline bool frozen();
            
            inline bool referred();
            
            
            
    // Member Operators


        //- Overridable function to handle the particle hitting a patch
        //  Executed before other patch-hitting functions
        bool hitPatch
        (
            const polyPatch&,
            trackingData& td,
            const label patchI,
            const scalar trackFraction,
            const tetIndices& tetIs
        );

        //- Overridable function to handle the particle hitting a processorPatch
        void hitProcessorPatch
        (
            const processorPolyPatch&,
            trackingData& td
        );

        //- Overridable function to handle the particle hitting a wallPatch
        void hitWallPatch
        (
            const wallPolyPatch&,
            trackingData& td,
            const tetIndices&
        );

        //- Overridable function to handle the particle hitting a polyPatch
        void hitPatch
        (
            const polyPatch&,
            trackingData& td
        );

        static void readFields(Cloud<agent>& mC);

        static void writeFields(const Cloud<agent>& mC);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const agent&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "agentI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
