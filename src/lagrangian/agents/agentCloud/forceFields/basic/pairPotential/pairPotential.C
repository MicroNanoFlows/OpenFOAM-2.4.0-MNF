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

#include "pairPotential.H"
#include "IFstream.H"
#include "graph.H"
#include "agentCloud.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pairPotential, 0);

defineRunTimeSelectionTable(pairPotential, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pairPotential::pairPotential
(
    agentCloud& cloud,
    const word& name,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(cloud.mesh())),
    cloud_(cloud),
    rU_(cloud_.redUnits()),
    name_(name),
    rCut_(readScalar(dict.lookup("rCut"))),
    rMin_(readScalar(dict.lookup("rMin"))),
    dr_(readScalar(dict.lookup("dr")))
{
    if(rU_.runReducedUnits())
    {
        rCut_ /= rU_.refLength();
        rMin_ /= rU_.refLength();
        dr_ /= rU_.refLength();
        rCutSqr_ = rCut_*rCut_;
    }
    
    // splitting the name using a delimeter "A-B" => "A" and "B"
    idList_.setSize(2);
    
//     Info << nl << "name = " << name_ << endl;
    
    std::string s = name_;
    std::string delimiter = "-";
    
    size_t pos = 0;
    std::string token;

    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        idList_[0]=token;
        s.erase(0, pos + delimiter.length());
        idList_[1]=s;
    }
    
/*    const cellZoneMesh& cellZones = mesh_.cellZones();
    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("pairPotential::pairPotential()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"controllersDict"
            << exit(FatalError);
    }
    
    if(dict.found("control"))
    {    
        control_ = Switch(dict.lookup("control"));
    }
//     readStateFromFile_ = Switch(dict.lookup("readStateFromFile"));*/
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<pairPotential> pairPotential::New
(
    agentCloud& cloud,
    const word& name,
    const dictionary& dict
)
{
    word pairPotentialName
    (
        dict.lookup("pairPotential")
    );

    Info<< "Selecting pairPotential "
         << pairPotentialName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(pairPotentialName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "pairPotential::New(const dictionary&) : " << endl
            << "    unknown pairPotential type "
            << pairPotentialName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<pairPotential>
	(
		cstrIter()(cloud, name, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pairPotential::~pairPotential()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// void pairPotential::updateStateControllerProperties
// (
//     const dictionary& newDict
// )
// {
// //     controllerDict_ = newDict.subDict("controllerProperties");
// 
//     //- you can reset the controlling zone from here. This essentially
//     //  means that the coupling zone can infact move arbitrarily. To make
//     //  this happen we probably need to devise a technique for automatically
//     //  changing the cellZone else where, and then calling this function to
//     //  reset the controlling zone in which the controller operates in.
// 
//     if (newDict.found("control"))
//     {
//         control_ = Switch(newDict.lookup("control"));
//     }
// /*
//     if (newDict.found("readStateFromFile"))
//     {
//         readStateFromFile_ = Switch(newDict.lookup("readStateFromFile"));
//     }*/
// }

// const labelList& pairPotential::controlZone() const
// {
//     return mesh_.cellZones()[regionId_];
// }
// 
// const word& pairPotential::regionName() const
// {
//     return regionName_;
// }
// 
// const bool& pairPotential::controlInterForces() const
// {
//     return controlInterForces_;
// }
// 
// bool& pairPotential::controlInterForces()
// {
//     return controlInterForces_;
// }
// 
// 
// const bool& pairPotential::writeInTimeDir() const
// {
//     return writeInTimeDir_;
// }
// 
// const bool& pairPotential::writeInCase() const
// {
//     return writeInCase_;
// }

const List<word>& pairPotential::idList() const
{
    return idList_;
}

void pairPotential::writeTables(const fileName& pathName)
{
    Info<< "Writing energy and force to file for potential "
        << name_ << endl;
            
    label nBins = label((rCut_ - rMin_)/dr_) + 1;            

    scalarField R(nBins, 0.0);
    scalarField U(nBins, 0.0);
    scalarField f(nBins, 0.0);
    
    for (label i=0; i<nBins; ++i)
    {
        scalar r = rMin_+dr_*i;
        R[i] = r;
        U[i] = energy(r);
        f[i] = force(r);
    }
    {
        OFstream file(pathName/name_+"-RU.xy");

        if(file.good())
        {
            forAll(U, i)
            {
                file 
                    << R[i] << "\t"
                    << U[i] << "\t"
                    << f[i]
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void pairPotentialModel::write()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    
    {
        OFstream file(pathName/name_+"-SI.xy");

        if(file.good())
        {
            forAll(U, i)
            {
                file 
                    << dr_*i*rU_.refLength() << "\t"
                    << U[i]*rU_.refEnergy() << "\t"
                    << f[i]*rU_.refForce()
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void pairPotentialModel::write()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }  
    } 
}



} // End namespace Foam

// ************************************************************************* //
