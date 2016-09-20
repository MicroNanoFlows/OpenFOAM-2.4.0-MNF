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

#include "borderModel.H"
#include "IFstream.H"
#include "graph.H"
#include "agentCloud.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(borderModel, 0);

defineRunTimeSelectionTable(borderModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
borderModel::borderModel
(
    Time& time,
    agentCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(cloud.mesh())),
    time_(time),
    cloud_(cloud),
    rU_(cloud_.redUnits())
//     name_(name),
//     rCut_(readScalar(dict.lookup("rCut"))),
//     rMin_(readScalar(dict.lookup("rMin"))),
//     dr_(readScalar(dict.lookup("dr")))
{
    
    borderList_ = List<vectorList>(dict.lookup("bordersList"));
    
    Info << "borderList = " << borderList_ << endl;
    
    // temporary
    treshold_ = 0.001;
    
    checkClosedEndedBorders();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<borderModel> borderModel::New
(
    Time& time,
    agentCloud& cloud,
    const dictionary& dict
)
{
    word borderModelName
    (
        dict.lookup("borderModel")
    );

    Info<< "Selecting borderModel "
         << borderModelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(borderModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "borderModel::New(const dictionary&) : " << endl
            << "    unknown borderModel type "
            << borderModelName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<borderModel>
	(
		cstrIter()(time, cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

borderModel::~borderModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void borderModel::checkClosedEndedBorders()
{
    forAll(borderList_, i)
    {
        const vector& rI = borderList_[i][0];
        const vector& rJ = borderList_[i][borderList_[i].size()-1];
        scalar rD = mag(rI - rJ);
        
        if(rD > treshold_)
        {
            FatalError
                << "borderModel::checkClosedEndedBorders() " << endl
                << "    This border has its end points not connected: " << nl
                << borderList_[i] << nl
                << "    Make the start point and end points the same." << endl;
        }
    }
}

const List<vectorList>& borderModel::borderList() const
{
    return borderList_;
}

void borderModel::writeBorders(const fileName& pathName)
{
 
}



} // End namespace Foam

// ************************************************************************* //
