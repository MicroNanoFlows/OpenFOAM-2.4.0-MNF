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

#include "coulombShifted.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(coulombShifted, 0);
addToRunTimeSelectionTable(pairPotentialModel, coulombShifted, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
coulombShifted::coulombShifted
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const reducedUnits& redUnits,
    const word& name, 
    const dictionary& dict
)
:
    pairPotentialModel(mesh, molCloud, redUnits, name, dict),
    constant_(1.0/(4.0 * constant::mathematical::pi * 8.854187817e-12))   
{
 
    if(redUnits.runReducedUnits())
    {
        constant_ = (1.0/(4.0 * constant::mathematical::pi * redUnits.epsilonPermittivity()));
    }
    else
    {
        constant_ = 1.0/(4.0*constant::mathematical::pi*8.854187817e-12);
    }
    
    useTables_ = false; 
    
    EB_ =  2/rCut_;
    EC_ = 1/(rCut_*rCut_);
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

coulombShifted::~coulombShifted()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar coulombShifted::unscaledEnergy(const scalar r) const
{
    return constant_/r;
}



scalar coulombShifted::force(const scalar r) const
{
    scalar force = constant_*( (1/(r*r)) - EC_);
    
    return force;
}
    
scalar coulombShifted::energy(const scalar r) const
{
    scalar energy = constant_*( (1/r) - EB_ + EC_*r);
    
    return energy;
}


const dictionary& coulombShifted::dict() const
{
    return pairPotentialProperties_;
}

void coulombShifted::write(const fileName& pathName)
{
//     Info<< "Writing energy and force to file for potential "
//             << name_ << endl;
//             
// //     label nBins = 100000;
//     label nBins = label((rCut_ - rMin_)/dr_) + 1;            
//     
//     scalarField U(nBins, 0.0);
//     scalarField f(nBins, 0.0);
//     
//     for (label i=0; i<nBins; ++i)
//     {
//         scalar r = rMin_+dr_*i;
//         
//         U[i] = energy(r);
//         f[i] = force(r);
//     }
//     {
//         OFstream file(pathName/name_+"-electrostatics-RU.xy");
// 
//         if(file.good())
//         {
//             forAll(U, i)
//             {
//                 file 
//                     << dr_*i << "\t"
//                     << U[i] << "\t"
//                     << f[i]
//                     << endl;
//             }
//         }
//         else
//         {
//             FatalErrorIn("void shortRangeElectrostatic::write()")
//                 << "Cannot open file " << file.name()
//                 << abort(FatalError);
//         }
//     }
//     
//     {
//         OFstream file(pathName/name_+"-electrostatics-SI.xy");
// 
//         if(file.good())
//         {
//             forAll(U, i)
//             {
//                 file 
//                     << dr_*i*rU_.refLength() << "\t"
//                     << U[i]*rU_.refEnergy() << "\t"
//                     << f[i]*rU_.refForce()
//                     << endl;
//             }
//         }
//         else
//         {
//             FatalErrorIn("void shortRangeElectrostatic::write()")
//                 << "Cannot open file " << file.name()
//                 << abort(FatalError);
//         }  
//     } 
}

} // End namespace Foam

// ************************************************************************* //
