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

#include "polyCNGraphenePeriodic.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"
#include "SortableList.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyCNGraphenePeriodic, 0);

addToRunTimeSelectionTable(polyConfiguration, polyCNGraphenePeriodic, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyCNGraphenePeriodic::polyCNGraphenePeriodic
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
//     const word& name
)
:
    polyConfiguration(molCloud, dict/*, name*/)
//     propsDict_(dict.subDict(typeName + "Properties"))
{
//     tNs_.clear();
    nCarbon_ = 0;
    nNitrogen_ = 0;
    nQ_ = 0;
    nP_ = 0;
    
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyCNGraphenePeriodic::~polyCNGraphenePeriodic()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polyCNGraphenePeriodic::setInitialConfiguration()
{
    scalar b = readScalar(mdInitialiseDict_.lookup("bondLength"));
    
    const reducedUnits& rU = molCloud_.redUnits();
    
    b /= rU.refLength();
    
    scalar theta = 120*constant::mathematical::pi/180.0;
    scalar c = 2.0*b*sin(0.5*theta);

    Info << nl << "Information on properties of graphene sheet" << nl << endl;
    
    Info << "c = " << c
         << nl << "b = " << b
         << endl;    

    boundBox box = mesh_.bounds();
    
    const scalar& bo = b;    
    
    vector nL = mdInitialiseDict_.lookup("lengthNormal");    
    vector nB = mdInitialiseDict_.lookup("breadthNormal");
    vector nZ = mdInitialiseDict_.lookup("heightNormal");    
    scalar Z = readScalar(mdInitialiseDict_.lookup("Z"));
    
    nL /= mag(nL);
    nB /= mag(nB);
    nZ /= mag(nZ);
    
    scalar Lx = box.span() & nL;
    scalar Ly = box.span() & nB;

    
    
    bool shift = false;

    vector Vshift = vector::zero;
    
    if(mdInitialiseDict_.found("shift"))
    {
        shift = Switch(mdInitialiseDict_.lookup("shift"));
        
        if(shift)
        {
            Vshift = b*nL;
        }
    }
    
    
//     label L = readLabel(mdInitialiseDict_.lookup("length"));
    scalar Lscal = Lx/(3.0*b);
    scalar Bscal = Ly/c;
    
    Info<< "L (scalar) = " << Lscal
        << ", B (scalar) = " << Bscal 
        << endl;

    label L = label(Lscal+0.5);
    label B = label(Bscal+0.5);
    
    scalar Xperiodic = L*3*b;
    scalar Yperiodic = B*c;    
    
    Info << "Set system to X = " << Xperiodic << " in dir = " << nL << endl;
    Info << "Set system to Y = " << Yperiodic << " in dir = " << nB << endl;
    
    scalar tol = 0.1;
    
    Info << "residual L = " << mag(L - Lscal) 
        << ", residual B = " << mag(B - Bscal) << endl;
    
    if(mag(L - Lscal) > tol)
    {    
        FatalErrorIn("polyCNGraphenePeriodicId::setInitialConfiguration()")
                << "Adjust mesh in x direction - residual = " << mag(L - Lscal)
                << exit(FatalError);
    }
    
    if(mag(B - Bscal) > tol)
    {    
        FatalErrorIn("polyCNGraphenePeriodicId::setInitialConfiguration()")
                << "Adjust mesh in y direction - residual = " << mag(B - Bscal)
                << exit(FatalError);
    }    

    
    // first layer of hexagon rings
    
    DynamicList<vector> Glayer;
    
    for (label i=1; i<L+1; i++)
    {
        vector r1 = (3*b*i)*nL;
        vector r2 = (b+(3*b*i))*nL;
        vector r3 = (1.5*b+(3*b*i))*nL + 0.5*c*nB;
        vector r4 = (b+(3*b*i))*nL + c*nB;
        vector r5 = (3*b*i)*nL + c*nB;
        vector r6 = ((3*b*i)-(0.5*b))*nL + 0.5*c*nB;
        
        Glayer.append(r1);
        Glayer.append(r2);
        Glayer.append(r3);
        Glayer.append(r4);
        Glayer.append(r5);
        Glayer.append(r6);
    }
    
    //Glayer.shrink();
    
    // shift entire sheet so that the first atom which was created 
    // becomes centred at (0.0, 0.0)
    
    forAll (Glayer, a)
    {
       Glayer[a] -= nL*3*b;
       Glayer[a] += nL*b;
    }    
   
    // make a template for the first layer of hexagons
    DynamicList<vector> Gsheet;
    
    forAll (Glayer, a)
    {
        Gsheet.append(Glayer[a]);
    }
    
    // make a copy of the layer in the other direction to create the sheet
    
    for (label j=1; j<B; j++)
    {
        forAll (Glayer, a)
        {
            Gsheet.append(Glayer[a] + c*j*nB);
        }
    }
    
    //Gsheet.shrink();
    
    // remove any atoms that overlap eachother (to prevent the MD code from blowing up)
    
    DynamicList<vector> Gnew;
    
    scalar tolerance = 0.1;
    
    forAll (Gsheet, i)
    {
        const vector& rI = Gsheet[i];
        
        bool overlapping = false;
        
        forAll (Gnew, j)
        {
            const vector& rJ = Gnew[j];
            scalar rMag = mag(rI - rJ);
        
            if (rMag < tolerance)
            {
                overlapping = true;
            }
        }
    
        if (!overlapping)
        {
            Gnew.append(rI);
        }
    }
    
    //Gnew.shrink();
    
    // shift sheet upwards so that there is a gap in the y direction
    
    forAll (Gnew, i)
    {
       Gnew[i] += nB*c/4;
    }  

    // shift sheet in the z direction
    
    forAll (Gnew, i)
    {
        Gnew[i] += nZ*Z;
    }
    
    // shift for stagerred approach
    if(shift)
    {
        forAll (Gnew, i)
        {
            Gnew[i] += Vshift;
            
            label cell = -1;
            label tetFace = -1;
            label tetPt = -1;

            mesh_.findCellFacePt
            (
                Gnew[i],
                cell,
                tetFace,
                tetPt
            );

            if(cell == -1)
            {
                Gnew[i] -= Lx*nL;
            }
        }        
    }
    
    // READ IN MORE PROPERTIES FROM DICTIONARY

    word molIdNName(mdInitialiseDict_.lookup("molIdNitrogen"));
    
    word molIdCName(mdInitialiseDict_.lookup("molIdCarbon"));

    const List<word>& idList(molCloud_.cP().molIds());

    label molIdN = findIndex(idList, molIdNName);

    if(molIdN == -1)
    {
        FatalErrorIn("polyCNGraphenePeriodicId::setInitialConfiguration()")
            << "Cannot find molecule id: " << molIdNName 
            << nl << "in moleculeProperties/idList."
            << exit(FatalError);
    }

    label molIdC = findIndex(idList, molIdCName);

    if(molIdC == -1)
    {
        FatalErrorIn("polyCNGraphenePeriodicId::setInitialConfiguration()")
            << "Cannot find molecule id: " << molIdCName 
            << nl << "in moleculeProperties/idList."
            << exit(FatalError);
    }

    scalar temperature = 0.0;
    vector bulkVelocity = vector::zero;

    if (mdInitialiseDict_.found("temperature"))
    {
        temperature = readScalar(mdInitialiseDict_.lookup("temperature"));
    }

    if (mdInitialiseDict_.found("bulkVelocity"))
    {
        bulkVelocity = mdInitialiseDict_.lookup("bulkVelocity");
    }
    
    initialiseVelocities_ = false;
    
    if(mdInitialiseDict_.found("initialiseVelocities"))
    {
        initialiseVelocities_ = Switch(mdInitialiseDict_.lookup("initialiseVelocities"));
    }
    
    
    

    bool tethered = false;

    //- assume graphene sheet molecules are frozen always
    bool frozen = true;

    if (mdInitialiseDict_.found("frozen"))
    {
        frozen = Switch(mdInitialiseDict_.lookup("frozen"));
    }
    
    label noCatomsCreated = 0;
    
    forAll(Gnew, i)
    {  
        vector p = Gnew[i];
        
        label cell = -1;
        label tetFace = -1;
        label tetPt = -1;

        mesh_.findCellFacePt
        (
            p,
            cell,
            tetFace,
            tetPt
        );

        if(cell != -1)
        {
            insertMolecule
            (
                p,
                cell,
                tetFace,
                tetPt,
                molIdC,
                tethered,
                frozen,
                temperature,
                bulkVelocity
            );

            noCatomsCreated++;
        }
        else
        {
            Info << "WARNING: molecule at position = " << p << ", out of mesh" << endl;
        }
    }


    Info << tab << " No of initial carbon added: " << noCatomsCreated << endl;        
    
    
    // Delete molecules at edge for periodicity

                
    DynamicList<polyMolecule*> edgeMoleculesRaw;
    
    {

        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        
        scalar rX = GREAT;
        
        for
        (
                mol = molCloud_.begin();
                mol != molCloud_.end();
                ++mol
        )
        {
            if(mol().id() == molIdC)
            {
                if(mol().position().x() < rX)
                {
                    rX = mol().position().x();
                }
            }
        }
        
        Info << "rX = " << rX << endl;
        
        for
        (
                mol = molCloud_.begin();
                mol != molCloud_.end();
                ++mol
        )
        {
            if(mol().id() == molIdC)
            {
                nCarbon_++;
                
                if(mol().position().x() < (rX + 0.05))
                {
                    edgeMoleculesRaw.append(&mol());
                }
            }
        }
    }
    
    Info << "no of edge molecules = " << edgeMoleculesRaw.size() << endl;
    
    // order edge molecules 
    
    DynamicList<polyMolecule*> edgeMolecules;
    
    label index = 0;
    edgeMolecules.append(edgeMoleculesRaw[index]);
    
    {
        DynamicList<polyMolecule*> edgeMoleculesTemp;
        
        forAll(edgeMoleculesRaw, i)
        {
            if(i != index)
            {
                edgeMoleculesTemp.append(edgeMoleculesRaw[i]);
            }
        }
        edgeMoleculesRaw.clear();
        edgeMoleculesRaw.transfer(edgeMoleculesTemp);
    }
    
    while(edgeMoleculesRaw.size() != 0)
    {
        vector rI = edgeMolecules[index]->position();
        scalar rD = GREAT;
        label idR = -1;

        forAll(edgeMoleculesRaw, i)
        {
            vector rJ = edgeMoleculesRaw[i]->position();
            
            scalar rIJMag = mag(rI- rJ);
            
            if(rIJMag <  rD)
            {
                idR = i;
                rD = rIJMag;
            }
        }
        
        edgeMolecules.append(edgeMoleculesRaw[idR]);
        
        {
            DynamicList<polyMolecule*> edgeMoleculesTemp;
            
            forAll(edgeMoleculesRaw, i)
            {
                if(i != idR)
                {
                    edgeMoleculesTemp.append(edgeMoleculesRaw[i]);
                }
            }
            edgeMoleculesRaw.clear();
            edgeMoleculesRaw.transfer(edgeMoleculesTemp);
        }
        
        index++;
    }
    
    // order molecules into a 2D matrix 
    
    List<DynamicList<polyMolecule*> > orderedList(edgeMolecules.size());

    Info << "starting molecules = "<< edgeMolecules.size() << endl;

    forAll(edgeMolecules, i)
    {
        Info << edgeMolecules[i]->position() << endl;
    }
    
    scalar tres=0.7;
    
    forAll(orderedList, i)
    {
        vector rI = edgeMolecules[i]->position();
        DynamicList<polyMolecule*> selectList;
        
        vector min = vector
                        (
                        rI.x()-1,
                        rI.y()-bo*tres,
                        rI.z()-bo*tres
                    );
                        
        vector max = vector
                        (
                        mesh_.bounds().max().x()+1,
                        rI.y()+bo*tres,
                        rI.z()+bo*tres                         
                        );
        
        boundedBox bb(min, max);
        
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        
        
        for
        (
                mol = molCloud_.begin();
                mol != molCloud_.end();
                ++mol
        )
        {
            if(mol().id() == molIdC)
            {
                if(bb.contains(mol().position()))
                {
                    polyMolecule* molI = &mol();
                    selectList.append(molI);
                }
            }
        }
        
        // order list
        
        SortableList<scalar> toSortList(selectList.size());
        
        forAll(toSortList, j)
        {
            toSortList[j]=selectList[j]->position().x();
        }
        
//         if (i == 0)
//         {
//             Info << "toSortList = " << toSortList << endl;
//         }
        
        toSortList.sort();
        
//             Info << " (after sort) toSortList = " << toSortList << endl;
        
        labelList indices = toSortList.indices();

//             Info << " indices = " << indices << endl;
        
        forAll(indices, j)
        {
            orderedList[i].append(selectList[indices[j]]);
        }
    }
    
    
    
    bool change1 = false;
    
    if (mdInitialiseDict_.found("change"))
    {
        change1 = Switch(mdInitialiseDict_.lookup("change"));
    }
    
    forAll(orderedList, i)
    {
        bool change2 = change1;
        
        forAll(orderedList[i], j)
        {
            if(change2)
            {
                nCarbon_--;
                nNitrogen_++;
                orderedList[i][j]->id() = molIdN;
                change2 = false;
            }
            else
            {
                change2 = true;
            }
        }
        
        if(change1)
        {
            change1 = false;
        }
        else
        {
            change1 = true;
        }
    }
    
        
            
    nQ_ = nNitrogen_;
            
    // number of molecules to add 
    DynamicList<polyMolecule*> molsToDel;        
    DynamicList<polyMolecule*> atomsInPyridinic; 
    
    {
        List<vector> pLocations = List<vector>(mdInitialiseDict_.lookup("pyridinicLocations"));    
        List<List<label> > choices = List<List<label> >(mdInitialiseDict_.lookup("choices"));             
        
        {
            forAll(pLocations, k)
            {
                label idI = -1;
                label idJ = -1;
                scalar rD = GREAT;
                
                // find the closest nitrogen atom
                
                forAll(orderedList, i)
                {
                    forAll(orderedList[i], j)
                    {
                        if (orderedList[i][j]->id() == molIdN)
                        {
                            vector rI = orderedList[i][j]->position();
                            
                            scalar rIJMAG = mag(rI - pLocations[k]);
                            
                            if(rIJMAG < rD)
                            {
                                rD = rIJMAG;
                                idI = i;
                                idJ = j;
                            }
                        }
                    }
                }
                /*
                // test which way we should go?
                
                vector vC1=neighbouringAtom(orderedList, idI, idJ+1, molIdC);
                vector vC2=neighbouringAtom(orderedList, idI, idJ-1, molIdC);
                vector vN = orderedList[idI][idJ]->position();
                vector c = (vC1+vC2+vN)/3.0;
                
                if( ((vC1-vN) & (c-vN)) > 0 )
                {*/
                    
                List<label> ringI(4);
                List<label> ringJ(4);
                
                ringI[0]=1;ringJ[0]=0;
                ringI[1]=2;ringJ[1]=0;
                ringI[2]=2;ringJ[2]=-1;
                ringI[3]=2;ringJ[3]=1;

                List<label> ringPI(18);
                List<label> ringPJ(18);
                
                ringPI[0]=0;ringPJ[0]=0;
                ringPI[1]=1;ringPJ[1]=-1;
                ringPI[2]=1;ringPJ[2]=1;
                ringPI[3]=2;ringPJ[3]=-2;
                ringPI[4]=2;ringPJ[4]=2;
                ringPI[5]=3;ringPJ[5]=-3;
                ringPI[6]=3;ringPJ[6]=-1;
                ringPI[7]=3;ringPJ[7]=1;
                ringPI[8]=3;ringPJ[8]=3;
                
                ringPI[9]=0;ringPJ[9]=-1;
                ringPI[10]=0;ringPJ[10]=1;
                ringPI[11]=1;ringPJ[11]=-2;
                ringPI[12]=1;ringPJ[12]=2;
                ringPI[13]=2;ringPJ[13]=-3;
                ringPI[14]=2;ringPJ[14]=3;
                ringPI[15]=3;ringPJ[15]=0;
                ringPI[16]=3;ringPJ[16]=-2;
                ringPI[17]=3;ringPJ[17]=2;
                
                
                if(findIndex(choices[k], 2) != -1)
                {
                    label dI = 0;
                    label dJ = 0;
                    
                    setPyridinic(molsToDel, atomsInPyridinic, orderedList, ringI, ringJ, ringPI, ringPJ, idI, idJ, dI, dJ);
                }
                if(findIndex(choices[k], 1) != -1)
                {
                    label dI = 0;
                    label dJ = -6;
                    
                    setPyridinic(molsToDel, atomsInPyridinic, orderedList, ringI, ringJ, ringPI, ringPJ, idI, idJ, dI, dJ);
                }                    
                if(findIndex(choices[k], 3) != -1)
                {
                    label dI = 0;
                    label dJ = 6;
                    
                    setPyridinic(molsToDel, atomsInPyridinic, orderedList, ringI, ringJ, ringPI, ringPJ, idI, idJ, dI, dJ);
                }
                if(findIndex(choices[k], 4) != -1)
                {
                    label dI = -3;
                    label dJ = -3;
                    
                    setPyridinic(molsToDel, atomsInPyridinic, orderedList, ringI, ringJ, ringPI, ringPJ, idI, idJ, dI, dJ);
                }
                if(findIndex(choices[k], 5) != -1)
                {
                    label dI = -3;
                    label dJ = 3;
                    
                    setPyridinic(molsToDel, atomsInPyridinic, orderedList, ringI, ringJ, ringPI, ringPJ, idI, idJ, dI, dJ);
                }
                if(findIndex(choices[k], 6) != -1)
                {
                    label dI = -6;
                    label dJ = 0;
                    
                    setPyridinic(molsToDel, atomsInPyridinic, orderedList, ringI, ringJ, ringPI, ringPJ, idI, idJ, dI, dJ);
                }                    
            }
        }
    }
    
    // get a ring of molecules surrounding the pyridinic rings to avoid switching those with carbon atoms
    
    DynamicList<polyMolecule*> temp;
    
    forAll(orderedList, i)
    {
        forAll(orderedList[i], j)
        {
            if(findIndex(atomsInPyridinic, orderedList[i][j]) != -1)
            {
                highlightAtomsAroundPyridinic(temp, atomsInPyridinic, orderedList, i, j+1);
                highlightAtomsAroundPyridinic(temp, atomsInPyridinic, orderedList, i, j-1);
                highlightAtomsAroundPyridinic(temp, atomsInPyridinic, orderedList, i+1, j);
                highlightAtomsAroundPyridinic(temp, atomsInPyridinic, orderedList, i-1, j);
            }
        }
    }
    
    forAll(temp, i)
    {
        atomsInPyridinic.append(temp[i]);
    }
    
    // bring back carbon atoms 
    label nCBack = readLabel(mdInitialiseDict_.lookup("noOfCarbonAtoms"));
    scalar prob = readScalar(mdInitialiseDict_.lookup("probability"));
    
    label count = 0;
    
    if(nCBack > 0)
    {
        while(count < nCBack)
        {
            IDLList<polyMolecule>::iterator mol(molCloud_.begin());
            
            for
            (
                    mol = molCloud_.begin();
                    mol != molCloud_.end();
                    ++mol
            )
            {
                if(count < nCBack)
                {
                    if(mol().id() == molIdN)
                    {
                        polyMolecule* molI = &mol();
                            
                        if(findIndex(atomsInPyridinic, molI) == -1)
                        {
                            if(molCloud_.rndGen().sample01<scalar>() <= prob)
                            {
                                molI->id()=molIdC;
                                
                                nCarbon_++;
                                nNitrogen_--;
                                count++;
                                nQ_--;
                            }
                        }
                    }
                }
            }        
        }
    }
//     Info << "count = "  << count << endl;


    Info<< " no of nitrogen atoms in pyridinic = " << atomsInPyridinic.size() << endl;
    
    
    forAll(molsToDel, m)
    {
        molCloud_.deleteParticle(*molsToDel[m]);
    }  
    
    scalar nC = scalar(nCarbon_);
    scalar nN = scalar(nNitrogen_);
    scalar nT = nC + nN;

    scalar nP = scalar(nP_);
    scalar nQ = scalar(nQ_);
    scalar nTN = nP + nQ;
    
    Info << "total number of atoms = " << nT << nl
            << "no. carbon atoms = " << nCarbon_ << "( " << (nC/nT)*100 << "% total)" << nl
            << "no. nitrogen atoms = " << nNitrogen_ << "( " << (nN/nT)*100 << "% total)" << nl
            << "no. P struct = " << nP  << "( " << (nP/nTN)*40 << "% total)" << nl
            << "no. Q struct = " << nQ  << "( " << (nQ/nTN)*40 << "% total)" << nl
            << endl;
             
      
}



void polyCNGraphenePeriodic::setPyridinic
(
    DynamicList<polyMolecule*>& molsToDel,    
    DynamicList<polyMolecule*>& atomsInPyridinic,
    const List<DynamicList<polyMolecule*> >& orderedList,
    const List<label>& ringI,
    const List<label>& ringJ,
    const List<label>& ringPI,
    const List<label>& ringPJ, 
    label idI,
    label idJ,
    label dI,
    label dJ
)
{
    forAll(ringI, i)
    {
        label I=idI+ringI[i]+dI;
        label J=idJ+ringJ[i]+dJ;
        removeMolecule(molsToDel, atomsInPyridinic, orderedList, I, J);
        highlightPyridinicAtoms(atomsInPyridinic, orderedList, I, J);
    }
    
    nCarbon_-= 3; nNitrogen_--; nP_+=6; nQ_-=7;
    
    forAll(ringPI, i)
    {
        label I=idI+ringPI[i]+dI;
        label J=idJ+ringPJ[i]+dJ;                            
        highlightPyridinicAtoms(atomsInPyridinic, orderedList, I, J);
    }
}

void polyCNGraphenePeriodic::highlightPyridinicAtoms
(
    DynamicList<polyMolecule*>& atomsInPyridinic,
    const List<DynamicList<polyMolecule*> >& orderedList,
    label idI,
    label idJ
)
{
    label i = idI;
    label j = idJ;
    
    label nSizeX=orderedList.size();
    
    if(i < 0)
    {
       i += nSizeX;
    }
    else if(i >= nSizeX)
    {
       i -= nSizeX;
    }
    
    label nSizeY = orderedList[i].size();
    
    if(j < 0)
    {
       j += nSizeY;
    }
    else if(j >= nSizeY)
    {
       j -= nSizeY;
    }       
    
//     orderedList[i][j]->fraction()=0.5;
    if(findIndex(atomsInPyridinic, orderedList[i][j]) == -1)
    {
        atomsInPyridinic.append(orderedList[i][j]);
    }
}

void polyCNGraphenePeriodic::highlightAtomsAroundPyridinic
(
    DynamicList<polyMolecule*>& atomsAroundPyridinic,
    DynamicList<polyMolecule*>& atomsInPyridinic,
    const List<DynamicList<polyMolecule*> >& orderedList,
    label idI,
    label idJ
)
{
    label i = idI;
    label j = idJ;
    
    label nSizeX=orderedList.size();
    
    if(i < 0)
    {
       i += nSizeX;
    }
    else if(i >= nSizeX)
    {
       i -= nSizeX;
    }
    
    label nSizeY = orderedList[i].size();
    
    if(j < 0)
    {
       j += nSizeY;
    }
    else if(j >= nSizeY)
    {
       j -= nSizeY;
    }       
    
//     orderedList[i][j]->fraction()=0.5;
    if(findIndex(atomsInPyridinic, orderedList[i][j]) == -1)
    {
        if(findIndex(atomsAroundPyridinic, orderedList[i][j]) == -1)
        {
            atomsAroundPyridinic.append(orderedList[i][j]);
        }
    }
}


void polyCNGraphenePeriodic::removeMolecule
(
    DynamicList<polyMolecule*>& molsToDel,
    DynamicList<polyMolecule*>& atomsInPyridinic,
    const List<DynamicList<polyMolecule*> >& orderedList,
    label idI,
    label idJ
)
{
    label i = idI;
    label j = idJ;
    
    label nSizeX=orderedList.size();
    
    if(i < 0)
    {
       i += nSizeX;
    }

    if(i >= nSizeX)
    {
       i -= nSizeX;
    }
    
    label nSizeY = orderedList[i].size();
    
    if(j < 0)
    {
       j += nSizeY;
    }

    if(j >= nSizeY)
    {
       j -= nSizeY;
    }    

//     orderedList[i][j]->fraction()=0.0;
    molsToDel.append(orderedList[i][j]);
    atomsInPyridinic.append(orderedList[i][j]);
}

void polyCNGraphenePeriodic::switchBackToCarbon
(
    const List<DynamicList<polyMolecule*> >& orderedList,
    label idI,
    label idJ,
    label molIdC
)
{
    label i = idI;
    label j = idJ;
    
    label nSizeX=orderedList.size();
    
    if(i < 0)
    {
       i += nSizeX;
    }

    if(i >= nSizeX)
    {
       i -= nSizeX;
    }
    
    label nSizeY = orderedList[i].size();
    
    if(j < 0)
    {
       j += nSizeY;
    }

    if(j >= nSizeY)
    {
       j -= nSizeY;
    }    

    
    if(orderedList[i][j]->id() != molIdC)
    {
        orderedList[i][j]->id() = molIdC;
        nCarbon_++;
        nNitrogen_--;
        nQ_--;
    }
}
vector polyCNGraphenePeriodic::neighbouringAtom
(
    const List<DynamicList<polyMolecule*> >& orderedList,
    label idI,
    label idJ,
    label molIdC
)
{
    label i = idI;
    label j = idJ;
    
    label nSizeX=orderedList.size();
    
    if(i < 0)
    {
       i += nSizeX;
    }

    if(i >= nSizeX)
    {
       i -= nSizeX;
    }
    
    label nSizeY = orderedList[i].size();
    
    if(j < 0)
    {
       j += nSizeY;
    }

    if(j >= nSizeY)
    {
       j -= nSizeY;
    }    

    vector v = vector::zero;
    
    if(orderedList[i][j]->id() == molIdC)
    {
        v = orderedList[i][j]->position();
    }
    else
    {
        Info << "ERROR" << endl;    
    }
    
    return v;
}



} // End namespace Foam

// ************************************************************************* //
