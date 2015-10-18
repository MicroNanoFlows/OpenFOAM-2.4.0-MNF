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

#include "pairPotentials.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// constructor
pairPotentials::pairPotentials
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const constantMoleculeProperties& cP,
    const reducedUnits& redUnits
)
:
    mesh_(mesh),
    cP_(cP),
    potentialsDict_
    (
        IOobject
        (
            "potentialDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    pairPotentialsList_(potentialsDict_.lookup("pairs")),
    pairNames_(pairPotentialsList_.size()),
    pairIds_(pairPotentialsList_.size()),
    pairPotentials_(pairPotentialsList_.size()),
    exclusions_(pairPotentialsList_.size())

{
   
    if(pairPotentialsList_.size() > 0 )
    {
        Info << nl << "Creating potentials: " << nl << endl;
    
        forAll(pairPotentialsList_, i)
        {
            const entry& potI = pairPotentialsList_[i];
            const dictionary& potIDict = potI.dict();
            const word& headerName=potI.keyword();
            
            pairPotentials_[i] = autoPtr<pairPotentialModel>
            (
                pairPotentialModel::New(mesh, molCloud, redUnits, headerName, potIDict)
            );
    
            pairNames_[i] = pairPotentials_[i]->type();
            pairIds_[i] = i;
            
            exclusions_[i] = pairPotentials_[i]->exclusions();
        }
        
        testPairPotentials();        
    }
    else
    {
         FatalErrorIn("pairPotentials.C") << nl
                << " Something went wrong. You have no potentials in your system/potentialDict"
                << nl << nl << "Check that the format is as follows: " 
                << " pairs " << nl
                << " ( " << nl
                << "    enter your pair potentials here "
                << " ); " 
                << nl << abort(FatalError);       
    }
    
   
    // electrostatic potential 
    
    word headerName = "electrostatic";
    
    dictionary& dict = potentialsDict_.subDict(headerName);
    
    electrostaticPotential_ = autoPtr<pairPotentialModel>
    (
        pairPotentialModel::New(mesh, molCloud, redUnits, headerName, dict)
    );


    
    rCut_ = maxRCut();
//     rCutSqr_ = rCut_*rCut_;
}


pairPotentials::~pairPotentials()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void pairPotentials::testPairPotentials()
{
    Info << nl << "Testing pair potentials..." << nl << endl;

    // test 1 - check that all pair sites have pair potential combinations
    
    List<word> pairSites = cP_.pairPotSiteIdList();
    
    Info << "cP_.pairPotSiteIdList() " << cP_.pairPotSiteIdList() << endl;
    
    // imp: first to set size of lists
    pairPotLinks_.setSize(pairSites.size());
    
    forAll(pairSites, i)
    {
        pairPotLinks_[i].setSize(pairSites.size(), -1); 
    }
    
    forAll(pairSites, i)
    {
        forAll(pairSites, j)
        {
            bool pairsFound = false;
            
            forAll(pairPotentials_, k)
            {
                const List<word>& ids = pairPotentials_[k]->idList();
                
                if
                (
                    ( (ids[0] == pairSites[i]) && (ids[1] == pairSites[j]) ) ||
                    ( (ids[1] == pairSites[i]) && (ids[0] == pairSites[j]) ) 
                )
                {
                    pairPotLinks_[i][j] = k;
                    pairPotLinks_[j][i] = k;
                    pairsFound = true;
                }
            }
            
            if(!pairsFound)
            {
                FatalErrorIn("pairPotentials.C") << nl
                        << " You did not specify the following pair combination => "
                        << pairSites[i] << "-" << pairSites[j]
                        << nl << " in the system/potentialDict."
                        << nl << abort(FatalError);                 
                
            }
        }
    }
}




scalar pairPotentials::maxRCut()
{
    scalar rCut = 0.0;
    
    forAll(pairPotentials_, i)
    {
        if(pairPotentials_[i]->rCut() > rCut)
        {
            rCut = pairPotentials_[i]->rCut();
        }
    }
    
    if(electrostaticPotential_->rCut() > rCut)
    {
        rCut = electrostaticPotential_->rCut();
    }    

    
    return rCut;
}        

scalar pairPotentials::force
(
    const label pairSiteA,
    const label pairSiteB,
    const scalar r
) const
{
    label k = pairPotLinks_[pairSiteA][pairSiteB];
    
    return pairPotentials_[k]->force(r);
}

scalar pairPotentials::energy
(
    const label pairSiteA,
    const label pairSiteB,
    const scalar r
) const
{
    label k = pairPotLinks_[pairSiteA][pairSiteB];
    
    return pairPotentials_[k]->energy(r);
}


bool pairPotentials::rCutSqr
(
    const label pairSiteA,
    const label pairSiteB,
    const scalar rIJSqr
) const
{
    label k = pairPotLinks_[pairSiteA][pairSiteB];
    
    if(rIJSqr <= pairPotentials_[k]->rCutSqr())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// bool pairPotentials::excludeSites
// (
//     const label pairSiteA,
//     const label pairSiteB
// ) const
// {
//     
//     label k = pairPotLinks_[pairSiteA][pairSiteB];
//     
//     if(!exclusions_[k])
//     {
//         return false;
//     }
//     else
//     {
//         return pairPotentials_[k]->exclModel().excludeSites();
//     }
// }

scalar pairPotentials::rMin
(
    const label pairSiteA,
    const label pairSiteB
) const
{
     label k = pairPotLinks_[pairSiteA][pairSiteB];
    
     return pairPotentials_[k]->rMin();
}

const List< List<label> >& pairPotentials::pairPotLinks() const
{
    return pairPotLinks_;
}

const label& pairPotentials::pairPotLink
(
    const label& pairSiteA,
    const label& pairSiteB    
) const
{
    return pairPotLinks_[pairSiteA][pairSiteB];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
