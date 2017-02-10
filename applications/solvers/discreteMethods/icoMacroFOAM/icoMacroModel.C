/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2008-2009 OpenCFD Ltd.
    \\/      M anipulation   |
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
    icoMacroModel

Description



\*---------------------------------------------------------------------------*/

#include "icoMacroModel.H"

namespace Foam
{

    
    //- Constructor 
icoMacroModel::icoMacroModel
(
    Time& t,
    const reducedUnits& rU,
    const dictionary& dict
)
:
    time_(t),
    rU_(rU),
//     dict(dict.subDict("macroProperties")),
    nMicro_(readLabel(dict.lookup("noOfMicroElements"))),
    N_(readLabel(dict.lookup("noOfMicroSteps"))),
    nIter_(readLabel(dict.lookup("nIter"))),
    iter_(0),
    L_(readScalar(dict.lookup("macroLength"))),
    Fs_(readScalar(dict.lookup("Fs"))),    
    rhoN_(readScalar(dict.lookup("rhoN"))),
    n_(dict.lookup("unitVector"))
    
{
    Info << nl << "initialising macro domain" << endl;  
    
    if( nMicro_ <= 0 )
    {
        FatalErrorIn("icoMacroModel::icoMacroModel()")
            << "Please define more than 1 micro sudbomain."
            << exit(FatalError);
    }
    
    // place check on odd micro subdomains
    deltaS_= L_/scalar(nMicro_);
    
    s_.setSize(nMicro_, 0.0);
    mDot_.setSize(nMicro_, 0.0);
    phi_.setSize(nMicro_, 0.0);
    f_.setSize(nMicro_, 0.0);
    forces_.setSize(nMicro_, vector::zero);
    k_.setSize(nMicro_, 0.0);
    
    forAll(s_, i)
    {
        s_[i] = i*deltaS_;
    }
    
    mDotMean_.setSize(nIter_, 0.0);
    iteration_.setSize(nIter_, 0.0);
    
    mDotAll_.setSize(nMicro_);
    phisAll_.setSize(nMicro_);
    phisCoeffsAll_.setSize(nMicro_);
    fAll_.setSize(nMicro_);
    
    forAll(phisAll_, i)
    {
        mDotAll_[i].setSize(nIter_,0.0);
        phisAll_[i].setSize(nIter_,0.0);
        phisCoeffsAll_[i].setSize(nIter_,0.0);
        fAll_[i].setSize(nIter_,0.0);
    }
    
    

    
    nPts_ = 10000;
    sSmooth_.setSize(nPts_);
    phisSmooth_.setSize(nPts_);
    scalar ds = L_/nPts_;
    
    forAll(phisSmooth_, k)
    {
        sSmooth_[k] = ds*k;
        phisSmooth_[k].setSize(nIter_,0.0);
    }
    
    fileName fieldPath(time_.path()/"multiscale");
    
    if (isDir(fieldPath))
    {
        rmDir(fieldPath);
    }    
    
    mkDir(fieldPath);    
    
    casePath_ = fieldPath;
    
    Info << "end" << endl;
}



icoMacroModel::~icoMacroModel()
{}


// solve the macroscpoic equations

void icoMacroModel::solve() 
{
    Info << nl << "iter:" << iter_ << " solving macro equation"
        << endl;


    if(iter_==0)
    {
        for(label i=0; i < nMicro_; i++)
        {
            k_[i] = mDot_[i]/f_[i];
        }
    }
    
    scalar pi = constant::mathematical::pi;
    
    // set matrix 
    
    simpleMatrix<scalar> luMatrix(nMicro_+1, 0.0, 0.0);
    
    label r = 0; // row counter
    
    // overall pressure drop integral equation
    
    luMatrix[r][0] = L_;
    
    for(label j=1; j <= (nMicro_-1)/2; j++)
    {
        luMatrix[r][(2*j)-1] = -(L_/(2*pi*j))*(1-Foam::cos(2*pi*j));
        luMatrix[r][(2*j)+1-1] = L_*(Foam::sin(2*pi*j))/(2*pi*j);
    }
    
    luMatrix.source()[r]=0.0;
    
    r++;
    
    for(label i=0; i < nMicro_; i++)
    {
        luMatrix[r][nMicro_]=1;
        luMatrix[r][0]=-k_[i]/rhoN_;
        
        for(label j=1; j <= (nMicro_-1)/2; j++)
        {
            luMatrix[r][(2*j)-1] = -(k_[i]/rhoN_)*Foam::sin(2*pi*s_[i]*j/L_);
            luMatrix[r][(2*j)+1-1] = -(k_[i]/rhoN_)*Foam::cos(2*pi*s_[i]*j/L_);
        }
        
        luMatrix.source()[r]=mDot_[i]-(k_[i]/rhoN_)*phi_[i];
        r++;
    }    
    
    //solve matrix
    scalarField matrixSoln = luMatrix.LUsolve(); 
    
    // new phi coefficients 
    List<scalar> phiCoeffs(nMicro_,0.0);
    List<scalar> phi(nMicro_,0.0);
    
    for(label i=0; i < nMicro_; i++)
    {
        phiCoeffs[i]=matrixSoln[i];
        phi[i] = matrixSoln[0];
        
        for(label j=1; j <= (nMicro_-1)/2; j++)
        {
            phi[i] += matrixSoln[(2*j)-1]*Foam::sin(2*pi*s_[i]*j/L_);
            phi[i] += matrixSoln[(2*j)+1-1]*Foam::cos(2*pi*s_[i]*j/L_);
        }
    }
    
    // measure stuff here for output
    
    for(label i=0; i < nMicro_; i++)
    {
        mDotAll_[i][iter_]=mDot_[i];
        fAll_[i][iter_]=f_[i];
        phisCoeffsAll_[i][iter_]=phiCoeffs[i];        
    }
    
    phi_ = phi; // important
    
    iter_++; // important
    
    mDotMean_[iter_] = matrixSoln[nMicro_];
    iteration_[iter_] = scalar(iter_);

    for(label i=0; i < nMicro_; i++)
    {
        phisAll_[i][iter_]=phi_[i];
    }

    // setting phi(s) smooth function
    
    for(label k=0; k < nPts_; k++)
    {
//         scalar& p = phisSmooth_[k][iter_];
        phisSmooth_[k][iter_] = matrixSoln[0];        
        
        for(label j=1; j <= (nMicro_-1)/2; j++)
        {
            phisSmooth_[k][iter_] += matrixSoln[(2*j)-1]*Foam::sin(2*pi*sSmooth_[k]*j/L_);
            phisSmooth_[k][iter_] += matrixSoln[(2*j)+1-1]*Foam::cos(2*pi*sSmooth_[k]*j/L_);
        }        
    }
}

List<vector>& icoMacroModel::setAndGetForces()
{
    // determine new forces
    for(label i=0; i < nMicro_; i++)
    {
        f_[i]=(Fs_ + phi_[i])/rhoN_;
        forces_[i] = f_[i]*n_;
    }
    
    return forces_;
}

void icoMacroModel::write()
{
    // mass flow rate
    
    scalarField mDotMean(1,mDotMean_[iter_]);
    scalarField iteration(1,iteration_[iter_]);
    
    // mDots
    
    writeTimeData
    (
        casePath_,
        "mean_mDot_vs_iteration.xy",
        iteration,
        mDotMean,
        true
    );    
    
    writeVariablesSideways
    (
        casePath_,
        "mDots_vs_iteration.xy",
        iteration_[iter_],
        mDot_
    );

    writeVariables
    (
        casePath_,
        "mDots_vs_s.xy",
        s_,
        mDotAll_
    );
    
    // phis 
    
    writeVariablesSideways
    (
        casePath_,
        "phis_vs_iteration.xy",
        iteration_[iter_],
        phi_
    );    
   
    writeVariables
    (
        casePath_,
        "phis_vs_s.xy",
        s_,
        phisAll_
    );
    
    writeVariables
    (
        casePath_,
        "phis_vs_s_smooth.xy",
        sSmooth_,
        phisSmooth_
    );    

    // forces
    writeVariablesSideways
    (
        casePath_,
        "forces_vs_iteration.xy",
        iteration_[iter_],
        f_
    );    
   
    writeVariables
    (
        casePath_,
        "forces_vs_s.xy",
        s_,
        fAll_
    ); 
    
    
    if(iter_ == 1)
    {
        // k
        
        writeTimeData
        (
            casePath_,
            "k.xy",
            s_,
            k_
        );    
        
        writeTimeData
        (
            casePath_,
            "s.xy",
            s_
        );         
    }
    
}


void icoMacroModel::writeVariables
(
    const fileName& pathName,
    const word& nameFile,    
    const scalarField& xData, // nMicro
    const List<scalarField>& yData // nMicro outer, iter inner
)
{
    OFstream file(pathName/nameFile);

    if(file.good())
    {
        forAll(xData, i)
        {
            file << xData[i] << "\t";
                
            forAll(yData[i], j)
            {
                file << yData[i][j] << "\t";
            }
            
            file << endl;
        }
    }
    else
    {
        FatalErrorIn("void writeTimeData::writeTimeData()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}


void icoMacroModel::writeVariablesSideways
(
    const fileName& pathName,
    const word& nameFile,    
    const scalar& xData,
    const scalarField& yData
)
{
    fileName fName(pathName/nameFile);

    std::ofstream file(fName.c_str(),ios_base::app);
    file.precision(11);

    if(file.is_open())
    {
        file << xData << " ";
        
        forAll(yData, n)
        {
            file << yData[n] << " ";
        }
        
        file << nl;
    }
    else
    {
        FatalErrorIn("void writeTimeData::writeTimeData()")
            << "Cannot open file " << fName
            << abort(FatalError);
    }

    file.close();
    
}


const label& icoMacroModel::nMicro() const
{
    return nMicro_;
}
        
const label& icoMacroModel::nIter() const
{
    return nIter_;
}

scalarField& icoMacroModel::massFlowRates()
{
    return mDot_;
}

// const List<vector>& icoMacroModel::forces() const
// {
//     return f_;
// }

// void icoMacroModel::pseudoMD() 
// {
//     Info << "running pseudo MD solver" << endl;
//     
//     scalar ds=4e-9;
//     scalar dz=4e-9;
//     scalar targetMassDensity=1000;
//     scalar molMass=2.99e-26;
//     
//     scalar Fs=0.4e-12*targetMassDensity/molMass;
//     scalar standardDev=0.1;
//     scalar mu=0.855e-3;
//     scalar Ls=5e-9;
//     
//     Random rndGen(clock::getTime());
//     
//     scalarField h(nMicro_, 0.0);
//     
//     h[0]=0.5e-8;
//     h[1]=0.362e-8;
//     h[2]=0.1382e-8;
//     h[3]=0.1382e-8;
//     h[4]=0.362e-8;
// 
//     scalarField A(nMicro_, 0.0);
//     
//     forAll(A, i)
//     {
//         A[i]=h[i]*dz;
//     }
//     
//     for(label i=0; i < nMicro_; i++)
//     {
//         scalarField u(nPts_, 0.0);
//         scalarField y(nPts_, 0.0);
// 
//         for(label k=0; k < nPts_; k++)
//         {
//             y[k]=(h[i]/nPts_)*k;
//             u[k]=-(1/(2*mu))*(Fs+phi_[i])*(y[k]*y[k]-h[i]*y[k] - h[i]*Ls)+standardDev*rndGen.scalar01();
//         }
//         
//         polynomialLeastSquaresFit p
//         (
//             y,
//             u, 
//             2
//         );
// 
//         scalarField coeffs= p.coeffs();
//         
//         scalar uMean= (1/A[i])*( coeffs[2]*h[i]*h[i]*h[i]*h[i]/4
//                         + coeffs[1]*h[i]*h[i]*h[i]/3
//                         + coeffs[0]*h[i]*h[i]/2 );
//         
//         mDot_[i] = targetMassDensity*uMean*A[i];
//     }
// }

}// End namespace Foam

