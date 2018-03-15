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

#include "polyReplaceMoleculesRandomly.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyReplaceMoleculesRandomly, 0);

addToRunTimeSelectionTable(polyConfiguration, polyReplaceMoleculesRandomly, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyReplaceMoleculesRandomly::polyReplaceMoleculesRandomly
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
//     const word& name
)
:
    polyConfiguration(molCloud, dict/*, name*/)
//     propsDict_(dict.subDict(typeName + "Properties"))
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyReplaceMoleculesRandomly::~polyReplaceMoleculesRandomly()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polyReplaceMoleculesRandomly::setInitialConfiguration()
{
    
    
    label initialSize = molCloud_.size();

    const scalar p(readScalar(mdInitialiseDict_.lookup("p")));
    
    //0 - 1
    

    
    const List<word>& idList(molCloud_.cP().molIds());

    word molIdName1(mdInitialiseDict_.lookup("molIdSelect"));
    
    word molIdName2(mdInitialiseDict_.lookup("molIdReplace"));
        
    label molId1 = findIndex(idList, molIdName1);

    if(molId1 == -1)
    {
        FatalErrorIn("polyReplaceMoleculesRandomly::setInitialConfiguration()")
            << "Cannot find molecule id 1: " << molIdName1 
            << nl << "in moleculeProperties/idList."
            << exit(FatalError);
    }

    label molId2 = findIndex(idList, molIdName2);

    if(molId2 == -1)
    {
        FatalErrorIn("polyReplaceMoleculesRandomly::setInitialConfiguration()")
            << "Cannot find molecule id 2: " << molIdName2 
            << nl << "in moleculeProperties/idList."
            << exit(FatalError);
    }
        
   
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
             
    label nAtoms = 0;
    for
    (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
    )
    {
        if(mol().id() == molId1)
        {
            nAtoms++;
        }
     
    }
    
    Info << "no of initial molecules = " << nAtoms << endl;

    //List<vector> positions_(nAtoms, vector::zero);
    //List<label> ids_(nAtoms, 0); 
    
    label nConvAtoms = 0;
    //label  c=0;
    
    for
    (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
    )
    {
        if(mol().id() == molId1)
        {
           scalar random = molCloud_.rndGen().sample01<scalar>();
           
           if(random <= p)
           {
                //positions_[c]=mol().position();
                mol().id() = molId2;
                nConvAtoms++;
           }
           
           //c++;
        }
     
    }
    
    Info << "converted atoms = " << nConvAtoms << endl;
    
    Info << "probability = " << nConvAtoms/nAtoms << endl; 
    
}



} // End namespace Foam

// ************************************************************************* //
