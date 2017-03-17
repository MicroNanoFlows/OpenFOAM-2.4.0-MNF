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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void forceFields::evaluatePair
(
    agent* molI,
    agent* molJ
)
{
  
    label idI = molI->id();
    label idJ = molJ->id();

    // pair potential interactions
//     controllers_.controlDuringForceComputation(molI, molJ); 

//     fields_.measurementsDuringForceComputation
//     (
//         molI,
//         molJ
//     );

    if(!molI->frozen() || !molJ->frozen())
    {
        // fraction
//         scalar f = molI->fraction();
//     
//         if(molJ->fraction() < f)
//         {
//             f = molJ->fraction();
//         }

        vector rIJ = molI->position() - molJ->position();
        scalar rIJMagSq = magSqr(rIJ);
        scalar rIJMag = mag(rIJ);
        
//         Info << "rIJMag = " << rIJMag  << endl;
        
        label k = pairPotentialLookUp_[idI][idJ];
                     
        if(rIJMagSq < pairPotentials_[k]->rCutSqr())
        {
            scalar energy = 0.0;
            vector force = vector::zero;

            pairPotentials_[k]->pairPotentialFunction
            (
                molI,
                molJ,
                rIJMag,
                energy,
                force
            );
            
//             vector fIJ = f*(rIJ/rIJMag)*force;

//             molI->f() += f*force;
//             molJ->f() += -f*force;

//             molI->potentialEnergy() += 0.5*f*energy;
//             molJ->potentialEnergy() += 0.5*f*energy;
            
            
//             Info<< molI->trackingNumber() << ", " << molJ->trackingNumber() 
//                 << ", rIJ = " << rIJMag << ", f = " << force << endl;
            
//             fields_.measurementsDuringForceComputationSite
//             (
//                 molI,
//                 molJ,
//                 sI,
//                 sJ
//             );
        }

        if(molI->R() > rIJMag)
        {
            molI->R() = rIJMag;
        }

        if(molJ->R() > rIJMag)
        {
            molJ->R() = rIJMag;
        }
    }
}


inline const scalar& forceFields::rCut() const
{
    return rCut_;
}

inline List< autoPtr<bodyForce> >& forceFields::bodyForces()
{
    return bodyForces_;
}

inline const List<vector>& forceFields::agentSocialForces() const
{
    return agentSocialForces_;
}

inline const List<vector>& forceFields::agentBodyForces() const
{
    return agentBodyForces_;
}
        
inline const List<vector>& forceFields::agentWallForces() const
{
    return agentWallForces_;
}

// inline List< autoPtr<agentController> >& agentControllers::stateControllers()
// {
// 	return stateControllers_;
// }
// 
// inline const List< autoPtr<agentController> >& agentControllers::stateControllers() const
// {
// 	return stateControllers_;
// }
// 
// inline const List<word>& agentControllers::stateControllersNames() const
// {
//     return sCNames_;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
