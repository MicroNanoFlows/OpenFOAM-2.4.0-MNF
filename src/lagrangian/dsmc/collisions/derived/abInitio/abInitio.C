/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "abInitio.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"
#include <sstream>

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(abInitio, 0);
    addToRunTimeSelectionTable(BinaryCollisionModel, abInitio, dictionary);
};



Foam::abInitio::abInitio
(
    const dictionary& dict,
    dsmcCloud& cloud
)
:
    BinaryCollisionModel(dict, cloud),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    xi_(),
    sigmaT_(),
    g_(),
    G_(readScalar(coeffDict_.lookup("G"))),
    nRows_(readLabel(coeffDict_.lookup("nRows"))),
    nColumns_(readLabel(coeffDict_.lookup("nColumns"))),
    xiTableFileName_(coeffDict_.lookup("deflectionAngleCosineTableFileName"))
{
    xi_.setSize(nRows_);
    
    forAll(xi_, i)
    {
        xi_[i].setSize(nColumns_ - 2);
    }
    
    sigmaT_.setSize(nRows_);
    g_.setSize(nRows_);
    
    //IFstream file(cloud.time().constant()/xiTableFileName_);
    
    word fileName = "./constant/" + xiTableFileName_;
    
    std::ifstream file(fileName.c_str());

    for(label row = 0; row < nRows_; row++)
    {
        std::string line;
        std::getline(file, line);
        if ( !file.good() )
            break;

        std::stringstream iss(line);

        for (label col = 0; col < nColumns_; col++)
        {
            std::string val;
            std::getline(iss, val, ',');

            std::stringstream convertor(val);
            
            if(col < nColumns_ - 2)
            {
                convertor >> xi_[row][col];
            }
            if(col == nColumns_ - 2)
            {
                convertor >> sigmaT_[row];
            }
            if(col == nColumns_ - 1)
            {
                convertor >> g_[row];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::abInitio::~abInitio()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::abInitio::active() const
{
    return true;
}



Foam::scalar Foam::abInitio::sigmaTcR
(
    const dsmcParcel& pP,
    const dsmcParcel& pQ
) const
{
    scalar cR = mag(pP.U() - pQ.U());

    label i = log(1.0 + cR/G_)/log(1.005);
    
    scalar sigmaTPQ = 0.0;
    
    if(i > 899)
    {
        i = 899;
    }
    
    if(i > 899 || i < 0)
    {
        FatalErrorIn("abInitio::sigmaTcR")
                << "Ab initio relative velocity out of scope."
                << exit(FatalError);
    }
    
    scalar g = g_[i];
    
    if(g > cR)
    {
        sigmaTPQ = sigmaT_[i];
        
        if(i != 0)
        {
            scalar x0 = g_[i-1];
            scalar x1 = g_[i];
            scalar y0 = sigmaT_[i-1];
            scalar y1 = sigmaT_[i];
            
            sigmaTPQ = y0 + (cR - x0)*((y1 - y0)/(x1 - x0));
        }
    }
    else if (g < cR)
    {
        Info << "g = " << g << endl;
        Info << "g+1 = " << g_[i+1] << endl;
        Info << "cR = " << cR << endl;
    }
    else
    {
        sigmaTPQ = sigmaT_[i];
    }

    //convert from square Angstroms  to square metres
    sigmaTPQ *= 1e-20; 

    return sigmaTPQ*cR;
}



void Foam::abInitio::collide
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    label& cellI
)
{
//     dsmcCloud& cloud_(this->owner());

    label typeIdP = pP.typeId();
    label typeIdQ = pQ.typeId();
    vector& UP = pP.U();
    vector& UQ = pQ.U();
    
    scalar collisionSeparation = sqrt(
            sqr(pP.position().x() - pQ.position().x()) +
            sqr(pP.position().y() - pQ.position().y())
    );
    
    cloud_.cellPropMeasurements().collisionSeparation()[cellI] += 
                                                        collisionSeparation;
    cloud_.cellPropMeasurements().nColls()[cellI]++;

    Random& rndGen(cloud_.rndGen());

    scalar mP = cloud_.constProps(typeIdP).mass();

    scalar mQ = cloud_.constProps(typeIdQ).mass();

    vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

    scalar cR = mag(UP - UQ);
    
    vector cRComponents = UP - UQ;

    label i = log(1.0 + cR/G_)/log(1.005);

    if(i > (nRows_ - 1))
    {
        i = (nRows_ - 1);
    }
    
    if(i > (nRows_ - 1) || i < 0)
    {
        FatalErrorIn("abInitio::collide")
                << "Ab initio relative velocity out of scope."
                << exit(FatalError);
    }
    
//     scalar g = g_[i];

    label j = rndGen.integer(0,99);
    
    if(j > 99 || j < 0)
    {
        FatalErrorIn("abInitio::collide")
                << "Ab initio cos(theta) out of scope."
                << exit(FatalError);
    }
    
    scalar cosTheta = xi_[i][j];
    
    if(i != 0)
    {
        scalar x0 = g_[i-1];
        scalar x1 = g_[i];
        scalar y0 = xi_[i-1][j];
        scalar y1 = xi_[i][j];
        
        cosTheta = y0 + (cR - x0)*((y1 - y0)/(x1 - x0));
    }
    

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*rndGen.scalar01();
    
    scalar D = sqrt(cRComponents.y()*cRComponents.y() + cRComponents.z()*cRComponents.z());
    
    scalar E = cR*cRComponents.z()*cos(phi);
    
    scalar F = cRComponents.x()*cRComponents.y()*sin(phi);
    
    scalar G = cR*cRComponents.y()*cos(phi);
    
    scalar H = cRComponents.x()*cRComponents.z()*sin(phi);

    vector postCollisionRelU =
       vector
        (
            cRComponents.x()*cosTheta + D*sin(phi)*sinTheta,
            cRComponents.y()*cosTheta + ((E - F)/D)*sinTheta,
            cRComponents.z()*cosTheta - ((G + H)/D)*sinTheta
        );


    UP = Ucm + postCollisionRelU*mQ/(mP + mQ);

    UQ = Ucm - postCollisionRelU*mP/(mP + mQ);
    
    label classificationP = pP.classification();
    label classificationQ = pQ.classification();
    
    //- class I molecule changes to class
    //- III molecule when it collides with either class II or class III
    //- molecules.
    
    if(classificationP == 0 && classificationQ == 1)
    {
        pP.classification() = 2;
    }
    
    if(classificationQ == 0 && classificationP == 1)
    {
        pQ.classification() = 2;
    }
    
    if(classificationP == 0 && classificationQ == 2)
    {
        pP.classification() = 2;
    }
    
    if(classificationQ == 0 && classificationP == 2)
    {
        pQ.classification() = 2;
    }
}

const Foam::dictionary& Foam::abInitio::coeffDict() const
{
    return coeffDict_;
}
// ************************************************************************* //
