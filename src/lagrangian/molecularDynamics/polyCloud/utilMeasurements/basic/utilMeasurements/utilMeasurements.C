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

#include "utilMeasurements.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

utilMeasurements::utilMeasurements
(
    Time& t,
    const polyMesh& mesh,
    const reducedUnits& rU    
)
:
    time_(t),
    dict_
    (
        IOobject
        (
            "inputDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    fieldList_(dict_.lookup("fields")),
    fieldNames_(fieldList_.size()),
    fieldIds_(fieldList_.size()),
    fields_(fieldList_.size())
{
    if(fields_.size() > 0)
    {
        fileName outputPath(time_.path()/"output");
        fileName inputPath(time_.path()/time_.timeName()/"uniform"/"poly");
       
        if (isDir(outputPath))
        {
            rmDir(outputPath);            
        }
        
        mkDir(outputPath);        

        forAll(fields_, f)
        {
            const entry& fieldI = fieldList_[f];
            const dictionary& fieldIDict = fieldI.dict();
    
            fields_[f] = autoPtr<utilField>
            (
                utilField::New(time_, mesh, rU, fieldIDict)
            );
    
            fieldNames_[f] = fields_[f]->type();
            fieldIds_[f] = f;

            fields_[f]->outputPath() = outputPath;
            fields_[f]->inputPath() = inputPath;
        }
    }
}


utilMeasurements::~utilMeasurements()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void utilMeasurements::createFields()
{
    forAll(fields_, f)
    {
        fields_[f]->createField();
    }
}

void utilMeasurements::calculateFields()
{
    forAll(fields_, f)
    {
        fields_[f]->calculateField();
    }
}

//- Note, not all fields automatically write out to disk.
void utilMeasurements::writeFields()
{
    fileName inputPath(time_.path()/time_.timeName()/"uniform"/"poly");
        
    forAll(fields_, f)
    {
        fields_[f]->inputPath() = inputPath;        
        fields_[f]->writeField();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
