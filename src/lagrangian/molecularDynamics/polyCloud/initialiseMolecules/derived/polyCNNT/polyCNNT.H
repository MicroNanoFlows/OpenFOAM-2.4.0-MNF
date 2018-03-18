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

Class
    polyCNNT

Description


SourceFiles
    polyCNNTI.H
    polyCNNT.C
    polyCNNTIO.C


\*---------------------------------------------------------------------------*/

#ifndef polyCNNT_H
#define polyCNNT_H
// 
#include "polyConfiguration.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyCNNT Declaration
\*---------------------------------------------------------------------------*/

class polyCNNT
:
	public polyConfiguration
{

private:

    // Private data

        //- subDictionary containing the properties
//         dictionary propsDict_;


    // Private Member Functions

        label nCarbon_;
        label nNitrogen_;

        label nQ_;
        label nP_;
        
        void removeMolecule
        (
            DynamicList<polyMolecule*>& molsToDel,
            DynamicList<polyMolecule*>& atomsInPyridinic,
            const List<DynamicList<polyMolecule*> >& orderedList,
            label idI,
            label idJ
        );
        
        void switchBackToCarbon
        (
            const List<DynamicList<polyMolecule*> >& orderedList,
            label idI,
            label idJ,
            label molIdC
        );
        
        void highlightPyridinicAtoms
        (
            DynamicList<polyMolecule*>& atomsInPyridinic,
            const List<DynamicList<polyMolecule*> >& orderedList,
            label idI,
            label idJ
        );

        vector neighbouringAtom
        (
            const List<DynamicList<polyMolecule*> >& orderedList,
            label idI,
            label idJ,
            label molIdC
        );         
        
public:

    //- Runtime type information
    TypeName("polyCNNT");

    // Constructors

        //- Construct from components
        polyCNNT
		(
            polyMoleculeCloud& molCloud,
            const dictionary& dict
//             const word& name
		);

    // Destructor

 		~polyCNNT();


    // Member Functions

        void setInitialConfiguration();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
