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

#include "agentProperties.H"
// #include "writeTimeData.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//-  Constructor
agentProperties::agentProperties
(
    const polyMesh& mesh
)
:
    mesh_(mesh)
{
    Info << nl << "Reading agentProperties dictionary." << endl;
    
    IOdictionary dict
    (
        IOobject
        (
            "agentProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    
    
    List<word> idList = dict.lookup("idList");    
    
    N_ = idList.size();
    
    names_.setSize(N_);
    
    velocityMax_.setSize(N_);
// 
//     desDir_.setSize(N_);
//     
//     desSpeed_.setSize(N_);    
    
//     masses_.setSize(N_);
    
    forAll(names_, i)
    {
        if(findIndex(names_, idList[i]) != -1)
        {
            FatalErrorIn("agentProperties.C") << nl
                    << " You have defined more than one time the agent id = "
                    << idList[i]
                    << " in agentProperties, idList = "
                    << idList
                    << ". All agent ids need to be unique"
                    << nl << abort(FatalError);            
        }
        
        names_[i] = idList[i];        
    }

    Info << "Info: number of agent type ids = " << N_ 
         << nl << "names = " << names_
         << endl;
         
    PtrList<entry> molList(dict.lookup("agentProperties"));

    if(molList.size() != N_)
    {
        FatalErrorIn("agentProperties.C") << nl
                << " Number of molecules in idList = " << N_
                << ", does not match number of agent entries"
                << " in agentProperties subdictionary = "
                << molList.size()
                << nl << abort(FatalError);
    }
    
    // reading in information from constant/moleculeProperties dictionary 
    forAll(molList, i)
    {
        const entry& molEntryI = molList[i];
        const dictionary& subDict = molEntryI.dict();
        
        if (molEntryI.keyword() != names_[i])
        {
            FatalErrorIn("constantMoleculeProperties.C")
                << "The header of the subdictionary = "  << molEntryI.keyword()
                << ", with id number = " << i 
                << ", does not match with the molecule id type = " << names_[i]
                << " in the full idList = "
                << names_
                << "WARNING: The order of the contents in the idList needs to "
                << " match the headers of the moleculeProperties subdictionaries."
                << nl << abort(FatalError);
        }        
        
        scalar velI = readScalar(subDict.lookup("maxVelocity"));
        velocityMax_[i]=velI;
        
        /*
        vector desDirI = subDict.lookup("desiredDirection");
        desDir_[i]=desDirI;        
        
        
        scalar desSpeedI = readScalar(subDict.lookup("desiredSpeed"));
        desSpeed_[i]=desSpeedI;  */     
        
    }
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentProperties::~agentProperties()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const label& agentProperties::nAgentTypes() const
{
    return N_;
}

const List<word>& agentProperties::agentIds() const
{
    return names_;
}

const List<scalar>& agentProperties::vMax() const
{
    return velocityMax_;
}
// 
// 
// const List<vector>& agentProperties::desDir() const
// {
//     return desDir_;
// }
// 
// const List<scalar>& agentProperties::desSpeed() const
// {
//     return desSpeed_;
// }


// const List<scalar>& agentProperties::mass() const
// {
//     return masses_;
// }    
// 
// const scalar& agentProperties::mass(const label& id) const
// {
//     return masses_[id];    
// }

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
