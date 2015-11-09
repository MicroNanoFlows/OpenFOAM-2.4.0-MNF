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

#include "FENE.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(FENE, 0);
addToRunTimeSelectionTable(pairPotentialModel, FENE, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
FENE::FENE
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud, 
    const reducedUnits& redUnits,
    const word& name, 
    const dictionary& dict
)
:
    pairPotentialModel(mesh, molCloud, redUnits, name, dict),
    propsDict_(dict.subDict(typeName + "Coeffs")),
    K_(readScalar(propsDict_.lookup("K"))),
    r0_(readScalar(propsDict_.lookup("r0")))      
{
    if(redUnits.runReducedUnits())
    {
        K_ /= (redUnits.refMass()/ (redUnits.refTime()*redUnits.refTime()));
        
        r0_ /= redUnits.refLength();
    }

    setLookupTables();    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

FENE::~FENE()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar FENE::unscaledEnergy(const scalar r) const
{
    scalar constant = 1 - (r/r0_)*(r/r0_);
    scalar exp = Foam::exp(exponent);
    
    return -0.5*K_*r0*r0*ln;
}


// bool FENE::read
// (
//     const dictionary& pairPotentialProperties,
//     const reducedUnits& rU
// )
// {
//     pairPotentialModel::read(pairPotentialProperties, rU);
// 
//     FENECoeffs_ = pairPotentialProperties.subDict(typeName + "Coeffs");
// 
//     FENECoeffs_.lookup("sigma") >> sigma_;
//     FENECoeffs_.lookup("epsilon") >> epsilon_;
// 
//     if(rU.runReducedUnits())
//     {
//         sigma_ /= rU.refLength();
//         epsilon_ /= rU.refEnergy();
//     }
// 
//     return true;
// }

const dictionary& FENE::dict() const
{
    return propsDict_;
}


} // End namespace Foam

// ************************************************************************* //
