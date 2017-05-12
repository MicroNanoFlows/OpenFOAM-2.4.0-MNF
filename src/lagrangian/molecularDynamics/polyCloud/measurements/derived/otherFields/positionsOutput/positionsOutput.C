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

#include "positionsOutput.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(positionsOutput, 0);
addToRunTimeSelectionTable(polyField, positionsOutput, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
positionsOutput::positionsOutput
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fields_(t, mesh, "dummy"),
    molIds_(),
    excludeSites_(),
    nameFile_(propsDict_.lookup("fileName"))
{
    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();
    

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

positionsOutput::~positionsOutput()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void positionsOutput::createField()
{
    selectSiteIds sites
    (
        molCloud_.cP(),
        propsDict_,
        "sitesToExclude"
    );

    List<word> siteNames = sites.siteIdNames();

    excludeSites_.transfer(siteNames);

    Info   << "sites to exclude: " << excludeSites_ << endl;
}

void positionsOutput::afterForce()
{}

void positionsOutput::calculateField()
{}

void positionsOutput::writeField()
{
    if(time_.outputTime())
    {
        fileName timePath(time_.path()/time_.timeName()/"uniform");
            
        DynamicList<vector> positions;
        
        {    
            IDLList<polyMolecule>::iterator mol(molCloud_.begin());

            for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
            {
                if(findIndex(molIds_, mol().id()) != -1)
                {
                    forAll(mol().sitePositions(), i)
                    {
                        if(findIndex(excludeSites_, molCloud_.cP().siteNames(mol().id())[i]) == -1)
                        {
                            positions.append(mol().position()); // error here
                        }
                    }
                }
            }
        }
        
        OFstream file(timePath/nameFile_);

        if(file.good())
        {
            forAll(positions, i)
            {
                file << positions[i].x() << " "
                    << positions[i].y() << " " 
                    << positions[i].z() << " " 
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void positionsOutput::writeField()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
}


void positionsOutput::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void positionsOutput::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& positionsOutput::fields() const
{
    return fields_;
}

} // End namespace Foam

// ************************************************************************* //
