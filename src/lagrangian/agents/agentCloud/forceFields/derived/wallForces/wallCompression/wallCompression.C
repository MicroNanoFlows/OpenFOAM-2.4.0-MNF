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
    wallCompression

Description

\*----------------------------------------------------------------------------*/

#include "wallCompression.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(wallCompression, 0);

addToRunTimeSelectionTable(agentWallForce, wallCompression, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
wallCompression::wallCompression
(
    Time& time,
    const dictionary& dict
)
:
    agentWallForce(time, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    A_(readScalar(propsDict_.lookup("A"))),
    kappa_(readScalar(propsDict_.lookup("kappa")))
{
    spaceVarying_ = true;
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

wallCompression::~wallCompression()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector wallCompression::force(agent* agentI, const vector& rW)
{
    scalar Ri = agentI->radius();                
    
    // vector between point-on-wall and agent's position
    vector rWI = agentI->position() - rW;
    
    scalar dWI = mag(rWI);
    
    vector force = vector::zero;
    
    if(dWI <= rCut_)
    {
        vector nij = rWI/dWI;
        vector tij = vector (-nij.y(), nij.x(), 0);
        
        if(dWI < 0.5*Ri)
        {
            force = 0.5*A_*(2-2*dWI) - kappa_*(Ri-dWI)*(agentI->v() & tij)*tij;
        }
        else
        {
            force = ( - 
                    kappa_*(Ri-dWI)*(agentI->v() & tij)*tij;
        }
    }
    
    return force;
}
            
vector wallCompression::force(agent* agentI, const vector& rW, const scalar& time)
{
    return vector::zero;
}
            
vector wallCompression::force(agent* agentI, const scalar& time)
{
    return vector::zero;
}

 
            
            
void wallCompression::updateForce()
{}
            
            
void wallCompression::write
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
