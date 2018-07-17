/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "allConfigurations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//- Null Constructor
allConfigurations::allConfigurations
(
//     const polyMesh& mesh,
    const IOdictionary& dict
)
:
    dict_(dict),
    configurationList_(),
    ids_(),
    configurations_()
{}


//- Constructor for mdInitialise
allConfigurations::allConfigurations
(
//     const polyMesh& mesh,
    molecules& cloud,
    const IOdictionary& dict    
)
:
    dict_(dict),
    configurationList_(dict_.lookup("models")),
    ids_(configurationList_.size()),
    configurations_(configurationList_.size())
{

    Info << nl << "Reading in all instructions" << nl << endl;
    
    if (configurations_.size() > 0)
    {
        forAll(configurations_, c)
        {
            const entry& configurationI = configurationList_[c];
            const dictionary& configurationIDict = configurationI.dict();

            configurations_[c] = autoPtr<configuration>
            (
                configuration::New(cloud, configurationIDict)
            );

            ids_[c] = c;
        }
    }
    else
    {
        Info << nl 
        << "OOPS. Something went wrong. It seems you didn't introduce any models." 
        << nl << " Try:" 
        << nl << nl 
        << "models" << nl
        << "(" << nl
        << "      model" << nl
        << "      { "<< nl
        << "          info from model1 goes here" << nl    
        << "      } "<< nl
        << "      model" << nl
        << "      { "<< nl
        << "          info from model2 goes here" << nl    
        << "      } "<< nl
        << "      model" << nl
        << "      { "<< nl
        << "          info from model3 goes here" << nl    
        << "      } "<< nl
        << ");" << nl << endl;
    }    
    
}


//- initial configuration

void allConfigurations::spawn()
{
    if(configurations_.size() > 0 )
    {
        Info << nl << "Let the spawning of molecules begin..." << nl << endl;
        
        forAll(configurations_, c)
        {
            Info << "Configuration #: " << c+1 << endl;
            
            configurations_[c]->spawn();
        }
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
