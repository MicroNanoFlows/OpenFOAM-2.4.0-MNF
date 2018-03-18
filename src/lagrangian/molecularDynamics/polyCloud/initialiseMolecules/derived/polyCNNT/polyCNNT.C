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

#include "polyCNNT.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyCNNT, 0);

addToRunTimeSelectionTable(polyConfiguration, polyCNNT, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyCNNT::polyCNNT
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
//     const word& name
)
:
    polyConfiguration(molCloud, dict/*, name*/)
//     propsDict_(dict.subDict(typeName + "Properties"))
{
    nCarbon_ = 0;
    nNitrogen_ = 0;
    nQ_ = 0;
    nP_ = 0;
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyCNNT::~polyCNNT()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polyCNNT::setInitialConfiguration()
{
    

    const reducedUnits& rU = molCloud_.redUnits();
    
    // READ IN PROPERTIES

    label N = readLabel(mdInitialiseDict_.lookup("N"));
    label M = readLabel(mdInitialiseDict_.lookup("M"));
    vector startPoint = mdInitialiseDict_.lookup("startPoint");
    vector endPoint = mdInitialiseDict_.lookup("endPoint");
    scalar bondLength = readScalar(mdInitialiseDict_.lookup("bondLengthSI"));
    
    bondLength /= rU.refLength();    
    
    vector cntNormal = (endPoint - startPoint) / mag(endPoint - startPoint);

    const scalar& bo = bondLength;

    scalar n = scalar(N);
    scalar m = scalar(M);

    word molIdNName(mdInitialiseDict_.lookup("molIdNitrogen"));
    
    word molIdCName(mdInitialiseDict_.lookup("molIdCarbon"));

    const List<word>& idList(molCloud_.cP().molIds());

    label molIdN = findIndex(idList, molIdNName);

    if(molIdN == -1)
    {
        FatalErrorIn("polyCNNTId::setInitialConfiguration()")
            << "Cannot find molecule id: " << molIdNName 
            << nl << "in moleculeProperties/idList."
            << exit(FatalError);
    }

    label molIdC = findIndex(idList, molIdCName);

    if(molIdC == -1)
    {
        FatalErrorIn("polyCNNTId::setInitialConfiguration()")
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

    bool tethered = false;

    bool stopBeforeCntRoll = false;

    //- ability to stop sheet before rolling, and view in its sheet format
    if (mdInitialiseDict_.found("stopCNTRoll"))
    {
        stopBeforeCntRoll = Switch(mdInitialiseDict_.lookup("stopCNTRoll"));
    }

    //- assume CNT molecules are frozen always
    bool frozen = true;

    if (mdInitialiseDict_.found("frozen"))
    {
        frozen = Switch(mdInitialiseDict_.lookup("frozen"));
    }

    // default value for triming the graphene sheet. Make sure it is a small number, 
    // smaller than the bond length

    scalar offset = 0.01;

    // over-right here to avoid compiling
    if (mdInitialiseDict_.found("trimSheetOffset"))
    {
        offset = readScalar(mdInitialiseDict_.lookup("trimSheetOffset"));
    }

    // you can specify a perpendicular vector 
    bool rotateAboutMidAxis = false;
    vector tPerp = vector::zero;

    if (mdInitialiseDict_.found("perpendicularVector"))
    {
        tPerp = mdInitialiseDict_.lookup("perpendicularVector");

        tPerp /= mag(tPerp);
        rotateAboutMidAxis = true;
        scalar rD = tPerp & cntNormal;

        if(rD > SMALL)
        {
            FatalErrorIn("void Foam::polyMoleculeCloud::createCNTs() : ")
                << "chosen perpendicularVector: " << tPerp 
                << " is not perpendicular to the cnt normal: " << cntNormal
                << ". Angle to normal is: " 
                << acos(rD)*180.0/constant::mathematical::pi 
                << " (degrees)."
                << nl
                << exit(FatalError);
        }
    }

    //Check that the chiral vector is correct
    if(N < M)
    {
        FatalErrorIn("void Foam::polyMoleculeCloud::createCNTs() : ")
            << " m is greater than n! "
            << exit(FatalError);
    }

    scalar length = mag(endPoint - startPoint);

    scalar theta = 
    atan
    (
        (sqrt(3.0)*n)/
        ((2.0*m)+ n)
    );

    vector a1 = vector(bo*sqrt(3.0), 0.0, 0.0);
    vector a2 = vector(0.5*bo*sqrt(3.0), 1.5*bo, 0.0);
    vector cH = n*a1 + m*a2;

//     scalar Ch = sqrt(3.0)*bondLength*sqrt(magSqr(n) + magSqr(m) + n*m);
//     scalar radius = (Ch*0.5)/constant::mathematical::pi;

    scalar D = mag(cH)/constant::mathematical::pi;

    // information 
    Info << nl << "Creating polyCNNT. " << nl << endl;

    Info<< "CNT INFORMATION: "<< nl
        << "Radius: " << D*0.5 << nl
        << "Diameter: " << D << nl
        << "Length: " << length << nl
        << "Theta: " << theta << nl
        << "Chiral vector, cH: " << cH << nl
        << "a1: " << a1 << nl
        << "a2: " << a2 << nl
        << endl;


    // Create standard graphene sheet:

    scalar hexX = sqrt(3.0)*bo;

    scalar LxLeft = length*sin(theta) + hexX;
    scalar LxRight = mag(cH)*cos(theta) + hexX;
    scalar Ly = length*cos(theta) + mag(cH)*sin(theta) + 2.0*bo;

    DynamicList<vector> firstLayerPositions(0);

    label nColsRight = LxRight/hexX;
    label nColsLeft = LxLeft/hexX;
    label nRows = Ly/hexX;

    for(label i=0; i<nColsRight; i++)
    {
        vector rA = vector(hexX*i, 0.0, 0.0);
        vector rB = vector(hexX*(0.5+i), 0.5*bo, 0.0);
        vector rC = vector(hexX*(0.5+i), 1.5*bo, 0.0);
        vector rD = vector(hexX*i, 2.0*bo, 0.0);

        firstLayerPositions.append(rA);
        firstLayerPositions.append(rB);
        firstLayerPositions.append(rC);
        firstLayerPositions.append(rD);
    }

    for(label i=1; i<nColsLeft; i++)
    {
        vector rA = vector(hexX*i, 0.0, 0.0) - vector(2.0*hexX*i, 0.0, 0.0);
        vector rB = vector(hexX*(0.5+i), 0.5*bo, 0.0) - vector(2.0*hexX*i, 0.0, 0.0);
        vector rC = vector(hexX*(0.5+i), 1.5*bo, 0.0) - vector(2.0*hexX*i, 0.0, 0.0);
        vector rD = vector(hexX*i, 2.0*bo, 0.0) - vector(2.0*hexX*i, 0.0, 0.0);

        firstLayerPositions.append(rA);
        firstLayerPositions.append(rB);
        firstLayerPositions.append(rC);
        firstLayerPositions.append(rD);
    }

    DynamicList<vector> grapheneSheetPositions(0);

    forAll(firstLayerPositions, p)
    {
        grapheneSheetPositions.append(firstLayerPositions(p));
    }

    //firstLayerPositions.shrink();
    
    for(label i=1; i<nRows; i++)
    {
        forAll(firstLayerPositions, p)
        {
            grapheneSheetPositions.append(firstLayerPositions[p] + vector(0.0, 3.0*bo*i, 0.0));
        }
    }

    //grapheneSheetPositions.shrink();

    // rotate graphene sheet

    const tensor rotationMatrix
    (
        cos(theta), -sin(theta), 0.0,
        sin(theta), cos(theta), 0.0,
        0.0, 0.0, 1.0
    );

    vectorField rotatedGrapheneSheet;

    rotatedGrapheneSheet.transfer(grapheneSheetPositions);

    forAll(rotatedGrapheneSheet, p)
    {
        rotatedGrapheneSheet[p] = transform(rotationMatrix.T(), rotatedGrapheneSheet[p]);
    }

    // truncate graphene sheet

    vector pMin = vector::zero;
    vector pMax = vector(mag(cH), length, 0.0);

    // adding  buffer (to avoid truncation errors)
    vector span = pMax - pMin;
    span /= mag(span);

    pMin -= span*offset;
    pMax += span*offset;

    boundBox bb(pMin, pMax);

    DynamicList<vector> truncatedGrapheneSheet(0);

    forAll(rotatedGrapheneSheet, p)
    {
        if(bb.contains(rotatedGrapheneSheet[p]))
        {
            truncatedGrapheneSheet.append(rotatedGrapheneSheet[p]);
        }
    }

    //truncatedGrapheneSheet.shrink();

    if(stopBeforeCntRoll)
    {
        Info << "WARNING: CNT NOT ROLLED!" << endl;

        label noCatomsCreated = 0;
        
        forAll(truncatedGrapheneSheet, c)
        {  
            vector p = truncatedGrapheneSheet[c];
            
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

            noCatomsCreated += 1;
        }
    
        if (Pstream::parRun())
        {
            reduce(noCatomsCreated, sumOp<label>());
        }
    
        Info << tab << "CNT molecules added: " << noCatomsCreated << endl;
    }
    else
    {
        // Roll cnt
        vectorField rolledSheet;
        
        rolledSheet.transfer(truncatedGrapheneSheet);

        scalar s = mag(cH); // circumference
        scalar R = 0.5*D; // radius

        forAll(rolledSheet, c)
        {
            scalar phi = rolledSheet[c].x()*2.0*constant::mathematical::pi/s;

            vector pDash = vector(R*sin(phi), rolledSheet[c].y(), R*cos(phi));

            rolledSheet[c] = pDash;
        }

        // remove overlaps

        DynamicList<label> overlappingMolecules(0);

        forAll(rolledSheet, c1)
        {
            const vector& p1 = rolledSheet[c1];

            forAll(rolledSheet, c2)
            {
                if(c2 > c1)
                {
                    const vector& p2 = rolledSheet[c2];
                    
                    scalar rD = mag(p2-p1);

                    if(rD < 0.5*bo)
                    {
                        if(findIndex(overlappingMolecules, c2) == -1)
                        {
                            overlappingMolecules.append(c2);
                        }
                    }
                }
            }
        }

        //overlappingMolecules.shrink();

        label Ncnt = (rolledSheet.size() - overlappingMolecules.size());

        vectorField cntMolecules(Ncnt, vector::zero);

        label counter = 0;

        forAll(rolledSheet, c)
        {
            if(findIndex(overlappingMolecules, c) == -1)
            {
                cntMolecules[counter] = rolledSheet[c];
                counter++;
            }
        }

        // orient CNT from the local co-ordinate axis used to setup CNT, to fit 
        // in the global co-ordinate axis where it was defined.


        vector t1 = vector::zero;
        vector t2 = vector::zero;

        if(rotateAboutMidAxis)
        {
            t1 = tPerp;
            t2 = cntNormal ^ t1;
        }
        else
        {
            scalar magV = 0.0;
            vector tangent;
        
            while (magV < SMALL)
            {
                vector testThis = molCloud_.rndGen().sampleVectorMD<vector>();
        
                tangent = testThis - (testThis & cntNormal)*cntNormal;
                magV = mag(tangent);
            }
    
            t1 = tangent/magV;
            t2 = cntNormal ^ t1;
        }

        Info << "Tangential vectors - t1: " << t1 << ", t2 : " << t2 << endl;

        forAll(cntMolecules, c)
        {
            cntMolecules[c] = startPoint + cntMolecules[c].y()*cntNormal 
                                + t1*cntMolecules[c].x() 
                                + t2*cntMolecules[c].z();
        }

        // create atoms

        label noCatomsCreated = 0;
        
        forAll(cntMolecules, c)
        {
            vector p = cntMolecules[c];

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
        
            if (cell != -1)
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

                noCatomsCreated += 1;
            }
            else
            {
                Info << "WARNING - Atom at specified position " << p
                     << ", does not correspond to a mesh cell.... deleting"
                     << nl
                     << endl;
            }
        }
        
        Info << tab << " No of initial carbon added: " << noCatomsCreated << endl;        
        
        
        // Delete molecules at edge for periodicity
        
                 
        DynamicList<polyMolecule*> edgeMoleculesRaw;
        
        {
            vector max = vector
                        (
                            mesh_.bounds().max().x()-0.2,
                            mesh_.bounds().max().y(),
                            mesh_.bounds().max().z()
                        );
            
            boundedBox bbMesh( mesh_.bounds().min(), max);
            
                    
            IDLList<polyMolecule>::iterator mol(molCloud_.begin());
            DynamicList<polyMolecule*> molsToDel;
            
            for
            (
                    mol = molCloud_.begin();
                    mol != molCloud_.end();
                    ++mol
            )
            {
                if(mol().id() == molIdC)
                {
                    if(!bbMesh.contains(mol().position()))
                    {
                        polyMolecule* molI = &mol();
                        molsToDel.append(molI);
                    }
                    else
                    {
                        nCarbon_++;
                        
                        if(mol().position().x() < 0.15)
                        {
                            edgeMoleculesRaw.append(&mol());
                        }
                    }
                }
            }
                        
            label nDeletedstage1 = molsToDel.size();
            
            forAll(molsToDel, m)
            {
                molCloud_.deleteParticle(*molsToDel[m]);
            }                    

            Info << "no of carbon atoms deleted due to periodicity "<< nDeletedstage1 << endl;
        }
        
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
        scalar tres=0.5;
        
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

//             Info << "toSortList = " << toSortList << endl;
            
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
             
        List<vector> pLocations = List<vector>(mdInitialiseDict_.lookup("pyridinicLocations"));    
//         label p = 3;
        
        DynamicList<polyMolecule*> molsToDel;        
        DynamicList<polyMolecule*> atomsInPyridinic; 
        
//         if(p==3)
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
                
                // test which way we should go?
                
                
                vector vC1=neighbouringAtom(orderedList, idI, idJ+1, molIdC);
                vector vC2=neighbouringAtom(orderedList, idI, idJ-1, molIdC);
                vector vN = orderedList[idI][idJ]->position();
                vector c = (vC1+vC2+vN)/3.0;
                
                if( ((vC1-vN) & (c-vN)) > 0 )
                {
                    
                    // rm ring 1 (top)
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI+1, idJ); 
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI+2, idJ); 
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI+2, idJ-1);
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI+2, idJ+1);

                    
                    
                    nCarbon_-= 3; nNitrogen_--; nP_+=6; nQ_-=7;
    
                    
                    // rm ring 2 (bottom left)
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI-1, idJ-2);
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI-1, idJ-3);
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI-1, idJ-4);
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI-2, idJ-3);
                
                    nCarbon_-= 3; nNitrogen_--; nP_+=6; nQ_-=7;                
                    
                    // rm ring 3 (bottom right)
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI-1, idJ+2);
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI-1, idJ+3);
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI-1, idJ+4);
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI-2, idJ+3);
                    
                    nCarbon_-= 3; nNitrogen_--; nP_+=6; nQ_-=7;                

                    
                    // central atom
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI, idJ);
                    
                    // bottom left
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI, idJ-2);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI, idJ-4);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI, idJ-6);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-1, idJ-1);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-1, idJ-5);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-2, idJ-2);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-2, idJ-4);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-3, idJ-3);
                    
                    //bottom right
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI, idJ+2);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI, idJ+4);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI, idJ+6);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-1, idJ+1);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-1, idJ+5);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-2, idJ+2);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-2, idJ+4);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-3, idJ+3);
                    
                    // top
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+1, idJ+1);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+1, idJ-1);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+2, idJ+2);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+2, idJ-2);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+3, idJ-3);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+3, idJ-1);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+3, idJ+1);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+3, idJ+3);
                    
                }
                else
                {
                    // rm ring 1 (bottom)
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI-1, idJ); 
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI-2, idJ); 
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI-2, idJ-1);
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI-2, idJ+1);
                    
                    nCarbon_-= 3; nNitrogen_--; nP_+=6; nQ_-=7;
    
                    
                    // rm ring 2 (top left)
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI+1, idJ-2);
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI+1, idJ-3);
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI+1, idJ-4);
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI+2, idJ-3);
                
                    nCarbon_-= 3; nNitrogen_--; nP_+=6; nQ_-=7;                
                    
                    // rm ring 3 (top right)
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI+1, idJ+2);
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI+1, idJ+3);
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI+1, idJ+4);
                    removeMolecule(molsToDel, atomsInPyridinic, orderedList, idI+2, idJ+3);
                    
                    nCarbon_-= 3; nNitrogen_--; nP_+=6; nQ_-=7;     

                    
                    // central atom
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI, idJ);
                    
                    // top left
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI, idJ-2);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI, idJ-4);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI, idJ-6);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+1, idJ-1);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+1, idJ-5);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+2, idJ-2);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+2, idJ-4);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+3, idJ-3);
                    
                    // top right
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI, idJ+2);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI, idJ+4);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI, idJ+6);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+1, idJ+1);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+1, idJ+5);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+2, idJ+2);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+2, idJ+4);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI+3, idJ+3);
                    
                    // bottom
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-1, idJ+1);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-1, idJ-1);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-2, idJ+2);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-2, idJ-2);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-3, idJ-3);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-3, idJ-1);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-3, idJ+1);
                    highlightPyridinicAtoms(atomsInPyridinic, orderedList, idI-3, idJ+3);                    
                }
            }
        }
        
        
        // bring back carbon atoms 
        label nCBack = readLabel(mdInitialiseDict_.lookup("noOfCarbonAtoms"));
        scalar prob = readScalar(mdInitialiseDict_.lookup("probability"));
        
        label count = 0;
        
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
        
        Info << "count = "  << count << endl;
        
        
            
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
}

void polyCNNT::highlightPyridinicAtoms
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
    
//     orderedList[i][j]->fraction()=0.5;
    atomsInPyridinic.append(orderedList[i][j]);
}

void polyCNNT::removeMolecule
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

void polyCNNT::switchBackToCarbon
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
vector polyCNNT::neighbouringAtom
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
