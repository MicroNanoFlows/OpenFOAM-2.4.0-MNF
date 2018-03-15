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

#include "polySwitchMoleculesRandomly.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polySwitchMoleculesRandomly, 0);

addToRunTimeSelectionTable(polyMolsToDeleteModel, polySwitchMoleculesRandomly, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polySwitchMoleculesRandomly::polySwitchMoleculesRandomly
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyMolsToDeleteModel(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    p_(readScalar(propsDict_.lookup("p")))
{
    
    {
        selectIds ids
        (
            molCloud_.cP(),
            propsDict_,
            "molIdSelect"
        );
        
        molId1_ = ids.molIds()[0];
    }
    
    {
        selectIds ids
        (
            molCloud_.cP(),
            propsDict_,
            "molIdReplace"
        );
        
        molId2_ = ids.molIds()[0];
    }    

    
    findMolsToDel();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polySwitchMoleculesRandomly::~polySwitchMoleculesRandomly()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void polySwitchMoleculesRandomly::findMolsToDel()
{
    label nAtoms = 0;
    
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
                

        for
        (
                mol = molCloud_.begin();
                mol != molCloud_.end();
                ++mol
        )
        {
            if(mol().id() == molId1_)
            {
                nAtoms++;
            }
        
        }
    }
    
    Info << "no of initial molecules = " << nAtoms << endl;

    //List<vector> positions_(nAtoms, vector::zero);
    //List<label> ids_(nAtoms, 0); 
    
    label nConvAtoms = 0;
    //label  c=0;

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    for
    (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
    )
    {
        if(mol().id() == molId1_)
        {
           scalar random = molCloud_.rndGen().sample01<scalar>();
           
           if(random <= p_)
           {
                //positions_[c]=mol().position();
                mol().id() = molId2_;
                nConvAtoms++;
           }

        }
     
    }
    
    Info << "converted atoms = " << nConvAtoms << endl;
    
    Info << "probability = " << scalar(nConvAtoms)/scalar(nAtoms) << endl; 
    
    // as a precaution: rebuild cell occupancy
    molCloud_.rebuildCellOccupancy();
    molCloud_.prepareInteractions();
}

} // End namespace Foam

// ************************************************************************* //
