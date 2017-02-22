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

Description

\*---------------------------------------------------------------------------*/

#include "agentConfiguration.H"
#include "IFstream.H"
#include "graph.H"
#include "agentCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(agentConfiguration, 0);
defineRunTimeSelectionTable(agentConfiguration, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
agentConfiguration::agentConfiguration
(
    agentCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(cloud.mesh())),
    cloud_(cloud),
    initialiseDict_(dict),
    nAgentsAdded_(0)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<agentConfiguration> agentConfiguration::New
(

    agentCloud& cloud,
    const dictionary& dict
)
{
    word agentConfigurationName
    (
        dict.lookup("model")
    );

    Info<< "Selecting agentConfiguration "
         << agentConfigurationName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(agentConfigurationName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "agentConfiguration::New(const dictionary&) : " << endl
            << "    unknown agentConfiguration type "
            << agentConfigurationName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<agentConfiguration>
	(
		cstrIter()(cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentConfiguration::~agentConfiguration()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// vector agentConfiguration::equipartitionLinearVelocity
// (
//     scalar temperature,
//     scalar mass
// )
// {
//     return sqrt(1.38064852e-23*temperature/mass)*vector
//     (
//         cloud_.rndGen().GaussNormal(),
//         cloud_.rndGen().GaussNormal(),
//         cloud_.rndGen().GaussNormal()
//     );
// }

vector agentConfiguration::setInitialVelocity
(
    const scalar& minV,
    const scalar& maxV 
)
{
    return vector
    (
        cloud_.rndGen().scalar01()*(maxV-minV) + minV,
        cloud_.rndGen().scalar01()*(maxV-minV) + minV,
        0.0
    );
}


scalar agentConfiguration::gaussianDistribution
(
    const scalar& mean,
    const scalar& range
 
)
{
    return mean + 0.25*range*cloud_.rndGen().GaussNormal();
}

// scalar agentConfiguration::uniformDistibution
// (
//     const scalar& mean,
//     const scalar& range
//  
// )
// {
//     return cloud_.rndGen().scalar01()*(maxV-minV) + minV;
// }



void agentConfiguration::insertAgent
(
    const point& position,
    const label cell,
    const label tetFace,
    const label tetPt, 
    const label& id,
    const scalar& mass, 
    const scalar& radius,
    const scalar& desiredSpeed,
//     const bool& tethered,
    const bool& frozen,
    const vector& v
)
{
    point specialPosition(vector::zero);

    label special = 0;

    if (frozen)
    {
        specialPosition = position;

        special = agent::SPECIAL_FROZEN;
    }

//     vector v = equipartitionLinearVelocity(temperature, cloud_.cP().mass(id));
// 
//     v += bulkVelocity;

           
    cloud_.createAgent
    (
        position,
        cell,
        tetFace,
        tetPt,     
        v,
        vector::zero,
        vector::zero,
        specialPosition,
        mass,
        radius,
        desiredSpeed,
        0.0,
        GREAT,
        1.0,
        special,
        id,
        0,
        cloud_.getTrackingNumber()
    );
}
/*
void agentConfiguration::insertMolecule
(
    const point& position,
    const label& id,
    const bool& tethered,
    const bool& frozen,
    const scalar& temperature,
    const vector& bulkVelocity
)
{
    label cell = -1;
    label tetFace = -1;
    label tetPt = -1;

    mesh_.findCellFacePt
    (
        position,
        cell,
        tetFace,
        tetPt
    );    
    
    if(cell != -1)
    {
        point specialPosition(vector::zero);

        label special = 0;

        if (tethered)
        {
            specialPosition = position;

            special = agent::SPECIAL_TETHERED;
        }

        if (frozen)
        {
            specialPosition = position;

            special = agent::SPECIAL_FROZEN;
        }

//         const agent::constantProperties& cP = cloud_.constProps(id);

        vector v = equipartitionLinearVelocity(temperature, cloud_.cP().mass(id));

        v += bulkVelocity;

        vector pi = vector::zero;

        tensor Q = I;

        if (!cloud_.cP().pointMolecule(id))
        {
            pi = equipartitionAngularMomentum(temperature, id);
            scalar phi(cloud_.rndGen().sample01<scalar>()*constant::mathematical::twoPi);
            scalar theta(cloud_.rndGen().sample01<scalar>()*constant::mathematical::twoPi);
            scalar psi(cloud_.rndGen().sample01<scalar>()*constant::mathematical::twoPi);

            Q = tensor
            (
                cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
                cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
                sin(psi)*sin(theta),
                - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
                - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
                cos(psi)*sin(theta),
                sin(theta)*sin(phi),
                - sin(theta)*cos(phi),
                cos(theta)
            );
        }

        cloud_.createMolecule
        (
            position,
            cell,
            tetFace,
            tetPt,     
            Q,
            v,
            vector::zero,
            pi,
            vector::zero,
            specialPosition,
            special,
            id,
            1.0,
            cloud_.getTrackingNumber()
        );
    }
    else
    {
        Info << "WARNING. Molecule not inserted since position out of mesh = "
            << position
            << endl;
    }
}*/

void agentConfiguration::deleteAgent
(
    agent& mol
)
{
    //- remove agent from cloud
    cloud_.deleteParticle(mol);
}


const label& agentConfiguration::nAgents() const
{
    return nAgentsAdded_;
}

} // End namespace Foam

// ************************************************************************* //
