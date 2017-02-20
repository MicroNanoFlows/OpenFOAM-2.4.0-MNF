/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    wallHelbingExponential

Description

\*----------------------------------------------------------------------------*/

#include "wallHelbingExponential.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(wallHelbingExponential, 0);

addToRunTimeSelectionTable(agentWallForce, wallHelbingExponential, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
wallHelbingExponential::wallHelbingExponential
(
    Time& time,
    const dictionary& dict
)
:
    agentWallForce(time, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
//     force_(propsDict_.lookup("force")),
    A_(readScalar(propsDict_.lookup("A"))),
    B_(readScalar(propsDict_.lookup("B"))),
    k_(readScalar(propsDict_.lookup("k"))),
    kappa_(readScalar(propsDict_.lookup("kappa"))),
    rCut_(readScalar(propsDict_.lookup("rCut")))
{
    spaceVarying_ = true;
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

wallHelbingExponential::~wallHelbingExponential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector wallHelbingExponential::force(agent* agentI, const vector& rW)
{
    scalar rI = agentI->radius();                
    
    // vector between point-on-wall and agent's position
    vector rWI = agentI->position() - rW;
    
    scalar dWI = mag(rWI);
    
    vector force = vector::zero;
    
    if(dWI <= rCut_)
    {
        vector nij = rWI/dWI;
        vector tij = vector (-nij.y(), nij.x(), 0);
        
        if(dWI >= rI)
        {
            force = A_*exp((rI-dWI)/B_)*nij;
            
//             Info << "rI = " << agentI->position()
//                 << ", rW = " << rW 
//                 << ", rWI = " << rW
//                 << ", dWI = " << dWI
//                 << ", force = " << force 
//                 << endl;
        }
        else
        {
            force = (A_*exp((rI-dWI)/B_) + k_*(rI-dWI))*nij - 
                    kappa_*(rI-dWI)*(agentI->v() & tij)*tij;
        }
    }
    
    return force;
}
            
vector wallHelbingExponential::force(agent* agentI, const vector& rW, const scalar& time)
{
    return vector::zero;
}
            
vector wallHelbingExponential::force(agent* agentI, const scalar& time)
{
    return vector::zero;
}

 
            
            
void wallHelbingExponential::updateForce()
{}
            
            
void wallHelbingExponential::write
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
