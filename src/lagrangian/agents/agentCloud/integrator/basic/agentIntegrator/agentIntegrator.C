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

#include "agentIntegrator.H"
#include "IFstream.H"
#include "graph.H"
#include "agentCloud.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(agentIntegrator, 0);

defineRunTimeSelectionTable(agentIntegrator, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
agentIntegrator::agentIntegrator
(
    Time& t,
    agentCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(cloud.mesh())),
    cloud_(cloud),
    time_(t)

{
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<agentIntegrator> agentIntegrator::New
(
    Time& t,
    agentCloud& cloud,
    const dictionary& dict
)
{
    word agentIntegratorName = "velocityVerlet";
    
    if(dict.found("integrator"))
    {
        const word agentIntegratorNameTemp
        (
            dict.lookup("integrator")
        );
        
        agentIntegratorName = agentIntegratorNameTemp;
    }

    Info<< "Selecting agentIntegrator "
         << agentIntegratorName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(agentIntegratorName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "agentIntegrator::New(const dictionary&) : " << endl
            << "    unknown agentIntegrator type "
            << agentIntegratorName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<agentIntegrator>
    (
        cstrIter()(t, cloud, dict)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentIntegrator::~agentIntegrator()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





} // End namespace Foam

// ************************************************************************* //
