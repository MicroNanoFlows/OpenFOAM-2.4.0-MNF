/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class
    molecules

Description

\*----------------------------------------------------------------------------*/

#include "molecules.H"


namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
molecules::molecules
(
//     const polyMesh& mesh,
    const IOdictionary& dict 
)
:
    dict_(dict)
{
    const List<word> names (dict.lookup("typeNames"));
    
    typeNames_.setSize(names.size());
    typeIds_.setSize(names.size());
    
    forAll(names, i)
    {
        typeNames_[i] = names[i];
        typeIds_[i] = i+1;
    }
    
    setBoundsBox(dict, bMesh_, "meshBoundsBox");   
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

molecules::~molecules()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void molecules::setBoundsBox
(
    const dictionary& propsDict,
    boundsBox& bb,
    const word& name 
)
{

    const dictionary& dict(propsDict.subDict(name));
    
    vector startPoint = dict.lookup("startPoint");
    vector endPoint = dict.lookup("endPoint");

    bb.resetBoundedBox(startPoint, endPoint);
    
    Info << "mesh vertices:" << endl;    
    
    Info << bb.min().x() << " " << bb.max().x() << " " << "xlo xhi" << endl;
    Info << bb.min().y() << " " << bb.max().y() << " " << "ylo yhi" << endl;
    Info << bb.min().z() << " " << bb.max().z() << " " << "zlo zhi"  <<endl;
}

void molecules::write()
{
    // writes .lammps trajectory
    Info << "Writing LAMMPS file" << endl;
    
//     Info << positions_ << endl;

    {
        OFstream file("lammps.traj");
        
        if(file.good())
        {
            file << positions_.size() << endl;
            file << endl;
            
            forAll(positions_, i)
            {
                file << types_[i] << " " 
                     << positions_[i].x() << " " 
                     << positions_[i].y() << " "
                     << positions_[i].z() << " " 
                     << endl;
            }
        }
    }
}



label molecules::size()
{
    return positions_.size();
}

DynamicList<vector>& molecules::positions()
{
    return positions_;
}

const DynamicList<vector>& molecules::positions() const
{
    return positions_;
}

DynamicList<label>& molecules::types()
{
    return types_;
}

const DynamicList<label>& molecules::types() const
{
    return types_;
}

DynamicList<scalar>& molecules::charges()
{
    return charges_;
}

const DynamicList<scalar>& molecules::charges() const
{
    return charges_;
}

const List<word>& molecules::typeNames() const
{
    return typeNames_;
}

const List<label>& molecules::typeIds() const
{
    return typeIds_;
}

const boundsBox& molecules::bMesh() const
{
    return bMesh_;
}

boundsBox& molecules::bMesh()
{
    return bMesh_;
}

label molecules::getType(const word& typeName)
{
    label id = findIndex(typeNames_, typeName);
    
    if(id != -1)
    {
        id = typeIds_[id];
    }
        
    return id;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
