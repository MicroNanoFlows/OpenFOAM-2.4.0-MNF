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

#include "dropletInfo.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "IFstream.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dropletInfo, 0);

addToRunTimeSelectionTable(utilField, dropletInfo, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dropletInfo::dropletInfo
(
    Time& t,
    const polyMesh& mesh,
    const reducedUnits& rU,
    const dictionary& dict
)
:
    utilField(t, mesh, rU, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.lookup("fieldName")),
    nameFile1_(propsDict_.lookup("nameFile_positions")),
    nameFile2_(propsDict_.lookup("nameFile_energies")),
    nameFile3_(propsDict_.lookup("nameFile_forces")),
    nameFile4_(propsDict_.lookup("nameFile_nPairs")),
    nameFile5_(propsDict_.lookup("nameFile_referenceMols")),
//     translate_(propsDict_.lookup("translate")),
//     rhoTolerance_(readScalar(propsDict_.lookup("rhoToleranceFrac"))),
    endTime_(readScalar(propsDict_.lookup("endTime"))),
    mass_(readScalar(propsDict_.lookup("mass")))
    
{
//     writeTime_ = GET;
    deltaT_ = time_.deltaT().value();
    
//     endTime_ = time_.endTime().value();
    
//     Info << "endTime = " << endTime_ << endl;
    
    
    mass_ /= rU_.refMass();
    
//     nSteps_ = endTime_/deltaT_;
    
//     nOutputSteps_=1000;
//     
//     if(propsDict_.found("nOutputSteps"))
//     {
//         nOutputSteps_ = readLabel(propsDict_.lookup("nOutputSteps"));
//     }
//     
//     nSegs_ = 5;
//     
//     if(propsDict_.found("nSegments"))
//     {
//         nSegs_ = readLabel(propsDict_.lookup("nSegments"));
//     }    
//     
    startTime_ = deltaT_*2.0;
    
    if(propsDict_.found("startTime"))    
    {        
        startTime_ = readScalar(propsDict_.lookup("startTime"));
    }
    
    dynamic_ = true;
    
    if(propsDict_.found("dynamic"))    
    {
        dynamic_ = Switch(propsDict_.lookup("dynamic")); 
    }
    
    nSteps_ = readLabel(propsDict_.lookup("nStepsOutput"));
    
    
//     val_ = 1;
    
    boundedBox bbMesh( mesh_.bounds().min(), mesh_.bounds().max() );
//     scalar Lx = bbMesh.span().x();
//     scalar Ly = bbMesh.span().x();
//     scalar Lz = bbMesh.span().z();
    
    box_ = bbMesh;
    
    lengthZ_ = box_.span().z();
    
/*    if(translate_ == "X")
    {
        translationVector_ = vector(1, 0, 0)*Lx;
    }
    
    if(translate_ == "Y")
    {
        translationVector_ = vector(0, 1, 0)*Ly;
    }       
    
    if(translate_ == "Z")
    {
        translationVector_ = vector(0, 0, 1)*Lz;
    }  */     

    
    initialiseMesh();

    initialiseBins();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dropletInfo::~dropletInfo()
{}





// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dropletInfo::createField()
{
    
}

void dropletInfo::reconstructDroplet()
{
    Info << " reconstruct droplet ... " << endl;
    
    {
        label N = 0;

        forAll(positions_, i)
        {
            if(fluidMols_[i])
            {
                N++;
            }
        }
        
        nFluidMols_ = N;

        Info << "no. of fluid mols = " << N << endl;

        List<vector> fluidPositions(N);
        List<label> fluidIds(N);
        
        label c = 0;
        
        forAll(positions_, i)
        {
            if(fluidMols_[i])
            {
                fluidPositions[c] = positions_[i];
                fluidIds[c] = i;
                c++;
            }
        }
        
        // populate mesh 
        
        forAll(fluidPositions, i)
        {
            const vector& rI = fluidPositions[i];
            
            if(box_.contains(rI))
            {
                scalar rDx = rI & unitVectorX_;
                scalar rDy = rI & unitVectorY_;
                label nX = label(rDx/binWidthMeshX_);
                label nY = label(rDy/binWidthMeshY_);
                
                if( (nX >= 0) && (nY >= 0) )
                {
                    if(nX == nBinsMeshX_) 
                    {
                        nX--;
                    }
            
                    if(nY == nBinsMeshY_) 
                    {
                        nY--;
                    }
                }
                
                rhoMesh_[nX][nY] += mass_;
            }
        }
        
//         check which cells are empty 
       
        forAll(cellOcc_, x)
        {
            forAll(cellOcc_[x], y)
            {
                rhoMesh_[x][y] /= volumeMeshCell_;
                
                cellOcc_[x][y] = 0;
                
                if(rhoMesh_[x][y] > rhoThreshold_)
                {
                   cellOcc_[x][y] = 1; 
                }
                
                //reset 
                rhoMesh_[x][y] = 0.0;
            }
        }
        
        dropletMolPositions_.clear();
        dropletMolIds_.clear();

        forAll(fluidPositions, i)
        {
            const vector& rI = fluidPositions[i];
            
            if(box_.contains(rI))
            {
                scalar rDx = rI & unitVectorX_;
                scalar rDy = rI & unitVectorY_;
                label nX = label(rDx/binWidthMeshX_);
                label nY = label(rDy/binWidthMeshY_);
                
                if( (nX >= 0) && (nY >= 0) )
                {
                    if(nX == nBinsMeshX_) 
                    {
                        nX--;
                    }
            
                    if(nY == nBinsMeshY_) 
                    {
                        nY--;
                    }
                }
                
                if(cellOcc_[nX][nY] > 0)
                {
                    dropletMolPositions_.append(rI);
                    dropletMolIds_.append(fluidIds[i]);
                }
            }
        }
        
        Info<< "number of molecules after sort = "
            << dropletMolPositions_.size() 
            << endl;
        
        
//         dropletMolPositions_.transfer(fluidPositions);
//         dropletMolIds_.transfer(fluidIds);
        
        
        
        
        /*
        
        
        List<DynamicList<vector> > selectPositions(nSegs_);
        List<DynamicList<label> > selectFluidIds(nSegs_);        
        List<DynamicList<label> > selectIds(nSegs_);

        // initial estimate
        
        List<label> molsChosen(N, -1);
        
        forAll(selectPositions, s)
        {
            vector com = vector::zero;
            
            if(s == 0)
            {
                com = fluidPositions[0];                
            }
            else
            {
                bool foundMol = false;
                
                forAll(molsChosen, i)
                {
                    if(molsChosen[i] == -1)
                    {
                        com = fluidPositions[i];
                        foundMol = true;
                        break;
                    }
                }
                
                if(!foundMol)
                {
                    break;
                }
            }
            
            scalar count = 1.0;
            
            label molsCount = 0;
            
            bool finish = false;
            
            while(!finish)
            {
                forAll(fluidPositions, j)
                {
                    const vector& rJ = fluidPositions[j];
                    
                    scalar rIJ = mag(com - rJ);
                
                    if(rIJ < count*maxR_)
                    {
                        selectPositions[s].append(rJ);
                        selectFluidIds[s].append(fluidIds[j]);
                        selectIds[s].append(j);
                    }
                }
                
                com = vector::zero;
                
                forAll(selectPositions[s], j)
                {
                    com += selectPositions[s][j];
                }
                
                label newMolsCount = selectPositions[s].size();
                
                com /= scalar(newMolsCount);
                
//                     Info << "com = " << com << endl;
            
                count += 1.0;
                
                if(newMolsCount == molsCount)
                {
                    finish = true;
                }
                else
                {
                    selectPositions[s].clear();
                    selectFluidIds[s].clear();
                    selectIds[s].clear();                    
                    molsCount = newMolsCount;
                }
            }
            
            // set remaining molecules
            forAll(selectIds[s], i)
            {
                molsChosen[selectIds[s][i]] = s;
            }
        }
        
        
        Info << "sorted lists = " << endl;
        
        forAll(selectPositions, s)
        {
            Info << selectPositions[s].size() << endl;
        }
        
        // reduce lists
        label noOfSegments = 0;
        
        forAll(selectPositions, s)
        {
            if(selectPositions[s].size() > val_)
            {
                noOfSegments++;
            }
        }
        
        List<DynamicList<vector> > sortedPositions(noOfSegments);
        List<DynamicList<label> > sortedIds(noOfSegments);                
            
        label cS = 0;
        
        forAll(selectPositions, s)
        {
            if(selectPositions[s].size() > val_)
            {
                sortedPositions[cS].transfer(selectPositions[s]);
                sortedIds[cS].transfer(selectFluidIds[s]);
                cS++;
            }
        }
        
        
        dropletMolPositions_.clear();
        dropletMolIds_.clear();

        if(noOfSegments == 1)
        {
            forAll(sortedPositions[0], i)
            {
                dropletMolPositions_.append(sortedPositions[0][i]);
                dropletMolIds_.append(sortedIds[0][i]);
            }             
        }
        else if(noOfSegments == 2)
        {
            vector com1 = vector::zero;
            
            forAll(sortedPositions[0], i)
            {
                com1 += sortedPositions[0][i];
            }

            com1 /= sortedPositions[0].size();
            
            vector com2 = vector::zero;
            
            forAll(sortedPositions[1], i)
            {
                com2 += sortedPositions[1][i];
            } 
            
            com2 /= sortedPositions[1].size();
            
            vector c21 = com1 - com2;
            
            scalar sgn = sign(translationVector_ & c21);

            forAll(sortedPositions[0], i)
            {
                dropletMolPositions_.append(sortedPositions[0][i]);
                dropletMolIds_.append(sortedIds[0][i]);
            }
            
            forAll(sortedPositions[1], i)
            {
                dropletMolPositions_.append(sortedPositions[1][i] + sgn*translationVector_);
                dropletMolIds_.append(sortedIds[1][i]);
            }                
        }
        else if(noOfSegments > 2)
        {
            // place code here for future reference
            
            FatalErrorIn("void dropletInfo::calculateField()")
                << "Did not code for more than 2 segements"
                << exit(FatalError);
        }*/
    }
}

void dropletInfo::calculateField()
{
    
    Info << " Calculate Field " << endl;
    
    if
    (
        (time_.timeOutputValue() > 0) && // skipping some known issues in first time step
        (time_.timeOutputValue() >= startTime_) &&
        (time_.timeOutputValue() <= endTime_)
    )
    {
        readFromStorage();
        
        reconstructDroplet();
        
        calculateBinnedProperties();
    }
}


// assumes a 2D bin for cylindrical droplets
void dropletInfo::calculateBinnedProperties()
{
    Info << "calculate 2D bin properties " << endl;
    
    vector com = vector::zero;
    
    forAll(dropletMolPositions_, i)
    {
        com += dropletMolPositions_[i];
    }
    
    com /= scalar(dropletMolPositions_.size());
    
    vector C1 = vector(com.x(), h_, 0.0);
    
    vector C2 = vector(com.x()+Lx_, h_+Ly_, box_.span().z()); 
    
    boundedBox box1(C1,C2);    
//     Info << "C = " << C << endl;
// 
//     Info << "C1 = " << C2 << endl;
    
    vector C3 = vector(com.x()-dX_, h_, 0.0);
    
    vector C4 = vector(com.x()+dX_, h_+Ly_, box_.span().z()); 
    
    boundedBox box2(C3,C4);

    vector C5 = vector(0, h_, 0.0);
    
    vector C6 = vector(Lx2_, Ly2_+h_, box_.span().z()); 
    
    boundedBox box3(C5,C6);

    
    forAll(dropletMolPositions_, i)
    {
        const vector& rI = dropletMolPositions_[i];
        
//         if(box1.contains(rI))
        {
            vector rSI = rI - C1;
            scalar rDx = mag(rSI & unitVectorX_);
            scalar rDy = mag(rSI & unitVectorY_);
            label nX = label(rDx/binWidthX_);
            label nY = label(rDy/binWidthY_);
            
            if( (nX >= 0) && (nY >= 0) )
            {
                if(nX == nBinsX_) 
                {
                    nX--;
                }
           
                if(nY == nBinsY_) 
                {
                    nY--;
                }
            }
            
            nMols_[nX][nY] += 1.0;
            pE_[nX][nY] += energies_[dropletMolIds_[i]];
        }
        
        if(box2.contains(rI)) // measurements of 1D profile along centre point of droplet
        {
            vector rSI = rI - C3;
            scalar rD = rSI & unitVectorY_;
            label n = label(rD/binWidthY1D_);
           
            if((n >= 0) && (n < nBinsY1D_))
            {
                nMolsY_[n] += 1.0;
                pEY_[n] += energies_[dropletMolIds_[i]];
            }
        }
        
        if(box3.contains(rI)) // dynamic droplet
        {
            vector rSI = rI - C5;
            scalar rDx = rSI & unitVectorX_;
            scalar rDy = rSI & unitVectorY_;
            label nX = label(rDx/binWidthX2_);
            label nY = label(rDy/binWidthY2_);
            
            if( (nX >= 0) && (nY >= 0) )
            {
                if(nX == nBinsX2_) 
                {
                    nX--;
                }
           
                if(nY == nBinsY2_) 
                {
                    nY--;
                }
            }
            
            nMols2_[nX][nY] += 1.0;

            nMols3_[nX][nY] += 1.0;
        }
    }
    
    averagingCounter_ += 1.0;
    
    averagingCounter2_ += 1.0;
    
    outputDynamicContourPlot();
}

void dropletInfo::outputDynamicContourPlot()
{
    if(averagingCounter2_ >= nSteps_)
    {
        if(dynamic_)
        {
            Info << "calculate dynamic plot " << endl;
            
            List<scalarField> rho(nBinsX2_);

            scalar volume = binWidthX2_*binWidthY2_*lengthZ_;
            
            forAll(rho, i)
            {
                rho[i].setSize(nBinsY2_, 0.0);
                
                forAll(rho[i], j)
                {
                    rho[i][j] = rU_.refMassDensity()*nMols2_[i][j]*mass_/(volume*averagingCounter2_);
                }
            }     
            
            outputCounter_++;
            std::string s;
            std::stringstream out;
            out << outputCounter_;
            s = out.str();        
            
            writeTimeData
            (
                outputPath_,
                s+"_density.txt",
                rho
            );            
            
            // reset 
            
            averagingCounter2_ = 0.0;
            
            forAll(nMols2_, i)
            {
                nMols2_[i] = 0.0;
            }
        }
    }
    
}

void dropletInfo::initialiseMesh()
{
    scalar Lx=box_.span().x();
    scalar Ly=box_.span().y();
    
    scalar binWidthMeshX = readScalar(propsDict_.lookup("binWidthMeshX"));
    scalar binWidthMeshY = readScalar(propsDict_.lookup("binWidthMeshY"));    
    
    //default
    
    binWidthMeshX_ = binWidthMeshX;
    binWidthMeshY_ = binWidthMeshY;
    
    nBinsMeshX_ = label(Lx/binWidthMeshX_) + 1;
    nBinsMeshY_ = label(Ly/binWidthMeshY_) + 1;
    
    binWidthMeshX_ = Lx/nBinsMeshX_;
    binWidthMeshY_ = Ly/nBinsMeshY_;
    
    Info << "MESH: binwidth X = " << binWidthMeshX_ 
         << ", binWidth Y = " << binWidthMeshY_ << endl;
    
    rhoMesh_.setSize(nBinsMeshX_);
    cellOcc_.setSize(nBinsMeshX_);
    
    forAll(rhoMesh_, x)
    {
        rhoMesh_[x].setSize(nBinsMeshY_, 0.0);
        cellOcc_[x].setSize(nBinsMeshY_, 0.0);
    } 
    
        
    volumeMeshCell_ = box_.span().z()*binWidthMeshX_*binWidthMeshY_;
        
    rhoThreshold_ = (1.0*mass_/volumeMeshCell_);    
    
    scalar NmolAv = (990/rU_.refMassDensity())*volumeMeshCell_/mass_;
    scalar NmolThreshold = rhoThreshold_*volumeMeshCell_/mass_;
    
    Info << "number of molecules expected per bin = " << NmolAv << endl;

    if(NmolAv <= NmolThreshold*10) 
    {
        FatalErrorIn("dropletInfo")
            << "it is recommended to increase size of cells "
            << " which are used to measure density "
            << " for triggering acceptance/rejection of mols for droplet."
            << exit(FatalError);        
    }
}

void dropletInfo::initialiseBins()
{
    unitVectorX_ = propsDict_.lookup("unitVectorX");
    unitVectorY_ = propsDict_.lookup("unitVectorY");
    unitVectorZ_ = propsDict_.lookup("unitVectorZ");
    
    nBinsX_ = readLabel(propsDict_.lookup("nBinsX"));
    nBinsY_ = readLabel(propsDict_.lookup("nBinsY"));
    h_ = readScalar(propsDict_.lookup("h"));
    Lx_ = readScalar(propsDict_.lookup("Lx"));    
    Ly_ = readScalar(propsDict_.lookup("Ly"));
    
    unitVectorX_ /= mag(unitVectorX_);
    unitVectorY_ /= mag(unitVectorY_);
    unitVectorZ_ /= mag(unitVectorZ_);  
    
    vector C1 = vector(0, h_, 0.0);
    
    vector C2 = vector(0+Lx_, h_+Ly_, box_.span().z());    
    
    boundedBox box(C1, C2);
    binWidthX_ = Lx_/nBinsX_;
    binWidthY_ = Ly_/nBinsY_;

        
    nMols_.setSize(nBinsX_);
    pE_.setSize(nBinsX_);
    
    forAll(nMols_, x)
    {
        nMols_[x].setSize(nBinsY_, 0.0);
        pE_[x].setSize(nBinsY_, 0.0);
    }    


    
    // 1D profile 
    dX_ = readScalar(propsDict_.lookup("dX"));
    nBinsY1D_ = readLabel(propsDict_.lookup("nBinsY1D"));
    binWidthY1D_ = Ly_/nBinsY1D_;
    
    nMolsY_.setSize(nBinsY1D_, 0.0);
    pEY_.setSize(nBinsY1D_, 0.0);
    
    
    Info << "binWidthY1D_ = " << binWidthY1D_ << endl;
    
    // dynamic
    
    nBinsX2_ = readLabel(propsDict_.lookup("nBinsX2"));
    nBinsY2_ = readLabel(propsDict_.lookup("nBinsY2"));    
    Lx2_ = readScalar(propsDict_.lookup("Lx2"));    
    Ly2_ = readScalar(propsDict_.lookup("Ly2"));

    vector C5 = vector(0.0, h_, 0.0);
    
    vector C6 = vector(Lx2_, Ly2_+h_, box_.span().z()); 
    
    boundedBox box3(C5,C6);
    binWidthX2_ = Lx2_/nBinsX2_;
    binWidthY2_ = Ly2_/nBinsY2_;    
    
    nMols2_.setSize(nBinsX2_);
    nMols3_.setSize(nBinsX2_);    
    forAll(nMols2_, x)
    {
        nMols2_[x].setSize(nBinsY2_, 0.0);
        nMols3_[x].setSize(nBinsY2_, 0.0);
    }    
    
    
    
    outputCounter_ = 0;
    averagingCounter_ = 0.0;
    averagingCounter2_ = 0;
}

    
void dropletInfo::writeField()
{
    // write xmol file
//     {
//         const reducedUnits& rU = rU_;
// 
//         fileName fName(outputPath_/time_.timeName()+"_shifted_"+fieldName_+".xmol");
//     
//         OFstream str(fName);
//     
//         str << dropletMolPositions_.size() << nl << "polyMoleculeCloud site positions in angstroms" << nl;
// 
//         // for all processors
//         forAll(dropletMolPositions_, i)
//         {
//             vector rI = dropletMolPositions_[i]*rU.refLength()*1e10;
// 
//             {
//                 str <<  'C'
//                     << ' ' << rI.x()
//                     << ' ' << rI.y()
//                     << ' ' << rI.z()
//                     << nl;
//             }
//         }
//     }    
//     
    if(time_.timeOutputValue() >= endTime_)
    {
        Info << "outputting... " << endl;
        
        List<scalarField> rho(nBinsX_);

        List<scalarField> bindingEnergy(nBinsX_);
        
        scalarField bindingEnergyY(nBinsY1D_, 0.0);
        
        forAll(pEY_, i)
        {
            if(nMolsY_[i] > 0.0)
            {
                bindingEnergyY[i] = rU_.refEnergy()*pEY_[i]/nMolsY_[i];
            }
        }
        
        scalar volume = binWidthX_*binWidthY_*lengthZ_;
        
        forAll(rho, i)
        {
            rho[i].setSize(nBinsY_, 0.0);
            bindingEnergy[i].setSize(nBinsY_, 0.0);
            
            forAll(rho[i], j)
            {
                rho[i][j] = rU_.refMassDensity()*nMols_[i][j]*mass_/(2.0*volume*averagingCounter_);
                
                if(nMols_[i][j] > 0.0)
                {    
                    bindingEnergy[i][j] = rU_.refEnergy()*pE_[i][j]/nMols_[i][j];
                }
            }
        }
        
        List<scalarField> rho3(nBinsX2_);
        scalar volume3 = binWidthX2_*binWidthY2_*lengthZ_;
            
        forAll(rho3, i)
        {
            rho3[i].setSize(nBinsY2_, 0.0);
            
            forAll(rho3[i], j)
            {
                rho3[i][j] = rU_.refMassDensity()*nMols3_[i][j]*mass_/(volume3*averagingCounter_);
            }
        } 
        
        scalarField binsX(nBinsX_, 0.0);
        
        forAll(binsX, i)
        {
            binsX[i] = (scalar(i) + 0.5)*binWidthX_*rU_.refLength();
        }
        
        scalarField binsY(nBinsY_, 0.0);
        
        forAll(binsY, i)
        {
            binsY[i] = (scalar(i) + 0.5)*binWidthY_*rU_.refLength();
        }
        
        scalarField binsY1D(nBinsY1D_, 0.0);
        
        forAll(binsY1D, i)
        {
            binsY1D[i] = (scalar(i) + 0.5)*binWidthY1D_*rU_.refLength();
        }       
        
        scalarField binsX2(nBinsX2_, 0.0);
        
        forAll(binsX2, i)
        {
            binsX2[i] = (scalar(i) + 0.5)*binWidthX2_*rU_.refLength();
        }
        
        scalarField binsY2(nBinsY2_, 0.0);
        
        forAll(binsY2, i)
        {
            binsY2[i] = (scalar(i) + 0.5)*binWidthY2_*rU_.refLength();
        }           
        

        writeTimeData
        (
            outputPath_,
            "binsX.txt",
            binsX
        );          


        writeTimeData
        (
            outputPath_,
            "binsY.txt",
            binsY
        ); 

        writeTimeData
        (
            outputPath_,
            "binsX2.txt",
            binsX2
        );          


        writeTimeData
        (
            outputPath_,
            "binsY2.txt",
            binsY2
        );         

        
        writeTimeData
        (
            outputPath_,
            "density_half.txt",
            rho
        );        

        writeTimeData
        (
            outputPath_,
            "density_full.txt",
            rho3
        );
        
        writeTimeData
        (
            outputPath_,
            "bindingEnergyY.txt",
            binsY1D,
            bindingEnergyY
        );
        
        writeTimeData
        (
            outputPath_,
            "bindingEnergy.txt",
            bindingEnergy
        );
        
        // write in gnuplot format
        
        {
            fileName fName(outputPath_/"density_gnuplot.txt");
        
            OFstream file(fName);
        
            file << "#x (A) \t y (A) \t rho g/cm^3" << nl;

            forAll(binsY, j)
            {
                scalar yI = (scalar(j) + 0.5)*binWidthY_*rU_.refLength()/1e-10;
                                
                forAll(binsX, i)
                {
                    scalar xI = (scalar(i) + 0.5)*binWidthX_*rU_.refLength()/1e-10;
                    
                    file << xI
                        << ' ' << yI
                        << ' ' << rho[i][j]*1000/(100*100*100)
                        << nl;                    
                    
                }                 
                
                file << nl;  
            }
        }
        
    }
}


void dropletInfo::readFromStorage()
{
    Info << "reading from storage" << endl;
    
    fileName pathName(inputPath_);
       
    {
        IFstream file(pathName/nameFile1_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> positions_;
        }
        else
        {
            FatalErrorIn("void dropletInfo::readFromStorage()")
                << "File: " << nameFile1_ << " not found"
                << exit(FatalError);            
        }
    }
    {
        IFstream file(pathName/nameFile2_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> energies_;
        }
        else
        {
            FatalErrorIn("void dropletInfo::readFromStorage()")
                << "File: " << nameFile2_ << " not found"
                << exit(FatalError);            
        }
        
    }    
    {
        IFstream file(pathName/nameFile3_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> forces_;
        }
        else
        {
            FatalErrorIn("void dropletInfo::readFromStorage()")
                << "File: " << nameFile3_ << " not found"
                << exit(FatalError);            
        }
        
    }
    {
        IFstream file(pathName/nameFile4_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> nPairs_;
        }
        else
        {
            FatalErrorIn("void dropletInfo::readFromStorage()")
                << "File: " << nameFile4_ << " not found"
                << exit(FatalError);            
        }
        
    }
    {
        IFstream file(pathName/nameFile5_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> fluidMols_;
        }
        else
        {
            FatalErrorIn("void dropletInfo::readFromStorage()")
                << "File: " << nameFile5_ << " not found"
                << exit(FatalError);            
        }
        
    }
}



} // End namespace Foam

// ************************************************************************* //
