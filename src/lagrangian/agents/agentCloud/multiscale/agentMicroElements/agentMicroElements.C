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

#include "agentMicroElements.H"
#include "simpleMatrix.H"
#include "polynomialLeastSquaresFit.H"
#include "fourierPolyLeastSquaresFit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    





//- Constructor 
agentMicroElements::agentMicroElements
(
    Time& t,
    const polyMesh& mesh,
    const agentProperties& cP
//     const dictionary& dict
)
:
    runTime_(t),
    cP_(cP),
    immDict_
    (
        IOobject
        (
            "immDict",
            runTime_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    
    propsDict_(immDict_.subDict("microElementProperties")),
    casePath_(t.path()/"measurements"),
    nMicro_(readLabel(propsDict_.lookup("noOfMicroElements"))),
    clouds_(nMicro_),
    N_(readLabel(propsDict_.lookup("noOfMicroSteps"))),
    macroLength_(readScalar(propsDict_.lookup("macroLength"))),
    macroDeltaT_(readScalar(propsDict_.lookup("macroDeltaT"))),     
    macroEndTime_(readScalar(propsDict_.lookup("macroEndTime"))), 
    macroTime_(t.startTime().value()),
    microTimes_(nMicro_, t.startTime().value()),
    microDeltaT_(t.deltaT().value()),
    deltaS_(0.0),
    s_(nMicro_, -1),
    nAgents_(nMicro_, 0.0),
    rhoM_(nMicro_, 0.0),
    m_(nMicro_, 0.0),
    W_(nMicro_, 0.0),

//     boxes_(),
//     bbVol_(nMicro_, 0.0),

    unitVector_(propsDict_.lookup("microUnitVector")),
    microLength_(readScalar(propsDict_.lookup("microLength"))),
    processorNumbers_(nMicro_, 0),
    procsMDaddressing_(Pstream::nProcs())
{
//     normal_ /= mag(normal_);
    
//     startFromZero_ = true;
//     
//     if(macroTime_ > 0.0)
//     {
//         startFromZero_ = false;
//     }
    
    Info << nl << " Creating sub-domains " << nl << endl;
    
    if( nMicro_ == 0 )
    {
        FatalErrorIn("agentMicroElements::agentMicroElements()")
            << "No molecule clouds defined. Define at least one cloud of molecules"
            << exit(FatalError);
    }
/*
    //- parallel-processing
    if(Pstream::parRun())
    {
        myProcNo_ = Pstream::myProcNo();
        
        parallelProcessingRequirements();
        
        //- creating chosen molecule clouds on current processor.
        if( molecules_.size() > 0 )
        {
            Pstream::removeParRun();

            forAll(molecules_, m)
            {
                //- to reduce memory allocation, only construct/initialise the assigned 
                // MD sub-domain(s) to the current processor
                if(Pstream::myProcNo() == processorNumbers_[m])
                {
                    std::string s;
                    std::stringstream out;
                    out << m;
                    s = out.str();            
                    
                    word cloudName = "atomisticMoleculeCloud_" + s;

                    Pout << nl << nl << "creating cloud #: " << m << ", name = " << cloudName << endl;
                    
                    molecules_[m] = autoPtr<atomisticMoleculeCloud>
                    (
                        atomisticMoleculeCloud::New(t, cloudName, mesh, pot, rU)
                    );
                }
            }

            Pstream::resetParRun();
        }        
    }*/

//     else // serial
//     {
//         if( clouds_.size() > 0 )
//         {
//             forAll(clouds_, m)
//             {
//                 std::string s;
//                 std::stringstream out;
//                 out << m;
//                 s = out.str();            
//                 
//                 word cloudName = "agentCloud_" + s;
// 
//                 Info << nl << nl << "creating cloud #: " << m << ", name = " << cloudName << endl;
//                 
//                 molecules_[m] = autoPtr<agentCloud>
//                 (
//                     agentCloud::New(t, mesh, cloudName, cP)
//                 );
//             }
//         }
//     }

    if( clouds_.size() > 0 )
    {
        forAll(clouds_, i)
        {
            std::string s;
            std::stringstream out;
            out << i;
            s = out.str();            
            
            word cloudName = "agentCloud_" + s;

            Info << nl << nl << "creating cloud #: " << i << ", name = " << cloudName << endl;
            
            clouds_[i] = autoPtr<agentCloud>
            (
                agentCloud::New(t, mesh, cloudName, cP, "NULL", true)
            );
            
//             initialise();
        }
    }
    
    // set macro time step
        
    Info << " macro time step = " << macroDeltaT_
         << ", micro time step = " << microDeltaT_
         << endl;
         
    unitVector_ /= mag(unitVector_);
    

    // set subdomain locations
    
    forAll(s_, i)
    {
        s_[i] = macroLength_*(i)/scalar(nMicro_-1);
    }

    deltaS_ = macroLength_/scalar(nMicro_-1);
    
    Info<< nl << " subdomain locations (s_i) = " << s_ 
        << nl << " subdomain spacing (delta s) = " << deltaS_
        << endl;

//     setBoundBoxes();
    
    // select liquid ids
    
//     liquidMolIds_.clear();
// 
//     selectIds ids
//     (
//         pot,
//         propsDict_,
//         "liquidMolIds"
//     );
// 
//     liquidMolIds_ = ids.molIds();  
    
    
    
//     Info << "areas = " << A_ << endl;
    
    
    // CAI //****
//     readVariableTimeSteppingInfo();
//     applyTimeStepping();
//     
//     setInsertionMethod();
    
//     readInEquationOfStates();
    
//     readInThermostatProperties();    
    
    
    // delete and remake directory for measuring output of properties
    
    fileName fieldPath(runTime_.path()/"measurements");

    if (!isDir(fieldPath))
    {
        mkDir(fieldPath);
    }

    // rule: clear fieldMeasurements every time the simulation is run

}


// void agentMicroElements::readVariableTimeSteppingInfo() 
// {
//     CAI_ = false;
//     
//     if (propsDict_.found("CAI"))
//     {
//         CAI_ = Switch(propsDict_.lookup("CAI"));   
//         
//         if(CAI_)
//         {
//             nMacro_ = (readLabel(propsDict_.lookup("nMacro")));
//             kG_ = (readScalar(propsDict_.lookup("kG")));
// 
//             Tmicro_ = (readScalar(propsDict_.lookup("Tmicro")));  
//             
//             {
//                 // read in period distribution
//                 
// //                 const word distributionName = propsDict_.lookup("periodName");
//                 const word distributionName = propsDict_.lookup("TmacroFileName");
// 
//                 fileName timePath(runTime_.system()/distributionName);
// 
//                 IFstream file(timePath);
// 
//                 List< Pair<scalar> > Tmacro;
// 
//                 if (file.good())
//                 {
//                     file >> Tmacro;
//                 }
//                 else
//                 {
//                     FatalErrorIn
//                     (
//                         "void agentMicroElements::readVariableTimeSteppingInfo()"
//                     )
//                         << "Cannot open file " << file.name()
//                         << abort(FatalError);
//                 }
// 
//                 label nBins = Tmacro.size();
//                 
//                 TMacro_.setSize(nBins, 0.0);
//                 timeMacro_.setSize(nBins, 0.0);
//                 
//                 forAll(Tmacro, bin)
//                 {
//                     timeMacro_[bin] = Tmacro[bin].first();                    
//                     TMacro_[bin] = Tmacro[bin].second();
//                 }
//                 
// //                 binWidth_ = period[1].first()-period[0].first();
//             }
//         }
//     }
// }

// void agentMicroElements::setTMacro() 
// {
//     // assume uniformly distributed points
//     scalar binWidth = timeMacro_[1] - timeMacro_[0];
//     
//     label index = label(macroTime_/binWidth);
//     
//     if(index < TMacro_.size()-1)
//     {
//         // linear interpolation
//         scalar y1 = TMacro_[index];
//         scalar y2 = TMacro_[index+1];
//         scalar x1 = timeMacro_[index];
//         scalar x2 = timeMacro_[index+1];
//         scalar x3 = macroTime_;
//         
//         Tmacro_ = y1 + ((x3-x1)*(y2-y1)/(x2-x1));
//     }
//     else
//     {
//         Info<< "WARNING in agentMicroElements::applyTimeStepping() " 
//             <<  nl << "exceeded the time-varying list. "
//             <<  endl;
//             
//         Tmacro_ = TMacro_[TMacro_.size()-1];
//     }
// }


/*
void agentMicroElements::applyTimeStepping() 
{
    if(CAI_)
    {
        //find the current Tmacro
        setTMacro(); 
        
        // find scale separation
        S_ = Tmacro_/Tmicro_;
        
        Info << "Scale separtion time-step = " << S_ << endl;
        
        // MODIFY MACRO TIME STEP
        macroDeltaT_ = Tmacro_/nMacro_;
        
        Info << "Macro time-step = " << macroDeltaT_ << endl;
        
        // MODIFY NO OF MICRO STEPS
       
        N_ = label( macroDeltaT_/(microDeltaT_*((kG_*(S_-1.0)) + 1.0)) );
        
        Info << "Micro steps (N) = " << N_ << endl;
        
        // g micro 
        gMicro_ = macroDeltaT_/(N_*microDeltaT_);
        
        Info << "Gmicro = " <<gMicro_ << endl;
    }
}*/
/*


void agentMicroElements::setNumericalIntegrationMethod()    
{
    const word schemeName = propsDict_.lookup("timeIntegration");
    
    if(schemeName == "euler")
    {
        euler_ = true;   
    }      
    else if(schemeName == "adamsBashforth")
    {
        adamsBashforth_ = true;             
    }
}*/

// void agentMicroElements::setNumericalDifferentiationMethod()    
// {
//     const word schemeName = propsDict_.lookup("differentiationScheme");
//     
//     if(periodic_)
//     {
//         if(schemeName == "three" )
//         {
//             three_ = true;   
//         }      
//         else if(schemeName == "five" )
//         {
//             five_ = true;   
//         }    
//         else if(schemeName == "seven" )
//         {
//             seven_ = true;   
//         }
//         else if(schemeName == "nine" )
//         {
//             nine_ = true;   
//         }
//         else if(schemeName == "fourier" )
//         {
//             fourier_ = true;
//             
//             // test if Nmicro are odd
//             if(nMicro_ % 2 == 0)
//             {
//                 // only fatal when nMicro is even 
//                 FatalErrorIn("agentMicroElements::agentMicroElements()")
//                     << "The fourier numerical differentiation only works for odd micro elements,"
//                     << " you specfied an even number = " 
//                     << nMicro_
//                     << exit(FatalError);               
//             }
//         }    
//         else if(schemeName == "fourierPolynomialSmoothing" )
//         {
//             fourierPolySmooth_ = true;
//             
//             label fourierOrder = readLabel(propsDict_.lookup("polynomialOrder"));
//             
//             // test if polynomialOrder_ is odd
//             if(fourierOrder % 2 == 0)
//             {            
//                 FatalErrorIn("agentMicroElements::agentMicroElements()")
//                     << "The fourier polynomial degree needs to be odd."
//                     << " You specfied an even number = " 
//                     << fourierOrder
//                     << exit(FatalError); 
//             }
//             
//             mFourierCoeffs_.setSize(fourierOrder, 0.0);
//             mOldFourierCoeffs_.setSize(fourierOrder, 0.0);
//             mOld2FourierCoeffs_.setSize(fourierOrder, 0.0);            
//             pFourierCoeffs_.setSize(nMicro_, 0.0);
//             rhoNFourierCoeffs_.setSize(nMicro_, 0.0);
//             
//             // test if Nmicro are odd
//             if(nMicro_ % 2 == 0)
//             {
//                 // only fatal when nMicro is even 
//                 FatalErrorIn("agentMicroElements::agentMicroElements()")
//                     << "The fourier numerical differentiation only works for odd micro elements,"
//                     << " you specfied an even number = " 
//                     << nMicro_
//                     << exit(FatalError);               
//             }            
//         }          
//         else 
//         {
//             FatalErrorIn("agentMicroElements::agentMicroElements()")
//                 << "No numerical differentiation scheme under the name = " 
//                 << schemeName
//                 << " for periodic techniques"
//                 << exit(FatalError);    
//         }        
//     }
//     else // non periodic
//     {
//         if(schemeName == "polynomial" )
//         {
//             polynomial_ = true;   
//         }
//         else if(schemeName == "polynomialSmoothing" )
//         {
//             polynomialSmoothing_ = true;   
//             polynomialOrder_ = readLabel(propsDict_.lookup("polynomialOrder"));
//             polynomialCoeffs_.setSize(polynomialOrder_+1);
//         }
//         else 
//         {
//             FatalErrorIn("agentMicroElements::agentMicroElements()")
//                 << "No numerical differentiation scheme under the name = " 
//                 << schemeName
//                 << " for non-periodic techniques"
//                 << exit(FatalError);    
//         }
//     }
// }


void agentMicroElements::updateMacroTime()
{
    macroTime_ += macroDeltaT_;
}

void agentMicroElements::updateMicroTime(const label& i)
{
    microTimes_[i] += microDeltaT_;
}

void agentMicroElements::updateMicroTimes()
{
    forAll(microTimes_, i)
    {
        microTimes_[i] += microDeltaT_;
    }
}

agentMicroElements::~agentMicroElements()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
/*
void agentMicroElements::getDeltaH(const scalar& rhoN, scalar& deltaH)
{
    label M = 1;
    
    scalarField coeffs(M+1,0.0);
    
//     coeffs[0] = -0.1136;
//     coeffs[1] = 0.6376;   
    
    coeffs[0] = -0.047203380767402;
    coeffs[1] = 0.346844932024265;
    
    deltaH = coeffs[0]*rhoN + coeffs[1];
}*/

/*
// not generally implemented
void agentMicroElements::updateHeights()
{
    forAll(boxes_, i)
    { 
        const boundBox& bb = boxes_[i];
        
        getDeltaH(rhoN_[i], deltaH_[i]);
        
        scalar hWall = bb.span().y();                
        
        h_[i] = hWall - (deltaH_[i]/frac_) + deltaH_[i];
        
        
//         if(nSteps_ > 0.0)
//         {
//             scalar hLiquid = liquidHeights_[i]/nSteps_;
//         
//             h_[i] = (bb.span().y() - hLiquid)*0.5 + hLiquid;            
//         }
        
        A_[i] = h_[i] * bb.span().z(); 
        
        bbVol_[i] = A_[i]*microLength_;
    }
    
    liquidHeights_ = 0.0;
    nSteps_ = 0.0;
    
//     Info << "heights" << endl;
}*/
/*
void agentMicroElements::thermostats()
{
    if(Pstream::parRun())
    {
        label p = Pstream::myProcNo();
        
        forAll(procsMDaddressing_[p], j)
        {
            label i = procsMDaddressing_[p][j];
            
            thermostat(i);
        }    
    }
    else
    {
        forAll(molecules_, i)
        {
            thermostat(i);
        }
    }
}*/
/*
void agentMicroElements::thermostat(const label& i)
{
    // measure mean velocity
    
    label nBins = tempBins_[i];
    
    scalarField mass(nBins, 0.0);
    vectorField mom(nBins, vector::zero);
    
    {
        IDLList<atomisticMolecule>::iterator mol(molecules_[i]->begin());

        for (mol = molecules_[i]->begin(); mol != molecules_[i]->end(); ++mol)
        {
            if(findIndex(liquidMolIds_, mol().id()) != -1)
            {
                const vector& rI = mol().position();
    
                label n = isPointWithinBin(rI, i);
    
                if(n != -1)
                {
                    const atomisticMolecule::constantProperties& constProp 
                                    = molecules_[i]->constProps(mol().id());
                                    
                    mass[n] += constProp.mass();
                    mom[n] += constProp.mass()*mol().v();
                }
            }
        }
    }
    
    vectorField velocity(nBins, vector::zero);
     
    forAll(mass, n)
    {
        if(mass[n] > 0.0)
        {
            velocity[n] = mom[n]/mass[n];
        }
    }

    // measure temperature
    
    scalarField kE(nBins, 0.0);
    scalarField dof(nBins, 0.0);
    
    {
        IDLList<atomisticMolecule>::iterator mol(molecules_[i]->begin());

        for (mol = molecules_[i]->begin(); mol != molecules_[i]->end(); ++mol)
        {
            if(findIndex(liquidMolIds_, mol().id()) != -1)
            {
                const vector& rI = mol().position();
    
                label n = isPointWithinBin(rI, i);
    
                if(n != -1)
                {
                    const atomisticMolecule::constantProperties& constProp 
                                    = molecules_[i]->constProps(mol().id());
                                    
                    kE[n] += 0.5*constProp.mass()*magSqr(mol().v() - velocity[n]);
                    dof[n] += constProp.degreesOfFreedom();
                }
            }
        }
    }

    scalarField measuredTemp(nBins, 0.0);
    scalarField chi(nBins, 1.0);     
    
    const scalar& kB = molecules_[i]->redUnits().kB();
    
    forAll(measuredTemp, n)
    {
        if(dof[n] > 0.0)
        {
            measuredTemp[n] = (2.0*kE[n])/(kB*dof[n]);

            chi[n] = sqrt(1.0 + (microDeltaT_/tauT_)*((temperature_/measuredTemp[n]) - 1.0) );
            
//             Info << "measured temperature = " <<  measuredTemp[n] << endl;
        }
    }
    
    // control velocities 
    
    {
        IDLList<atomisticMolecule>::iterator mol(molecules_[i]->begin());

        for (mol = molecules_[i]->begin(); mol != molecules_[i]->end(); ++mol)
        {
            if(findIndex(liquidMolIds_, mol().id()) != -1)
            {
                const vector& rI = mol().position();
    
                label n = isPointWithinBin(rI, i);
    
                if(n != -1)
                {
                    mol().v() = velocity[n] + (chi[n]*(mol().v() - velocity[n]));
                }
            }
        }
    }   
}*/
/*
void agentMicroElements::applyFourierSmoothingToMassFlux()
{
    // time smoothing
    
    if(PODmethod1_)
    {
        scalarField newS(nMicro_+1, 0.0);
        scalarField newM(nMicro_+1, 0.0);
        
        forAll(s_, i)
        {
            newS[i] = s_[i];
            newM[i] = m_[i];
        }
    
        newS[nMicro_]=macroLength_;        
        newM[nMicro_]=newM[0]; 
        
        label noSnapShots = noSnapShots_;
        
        Info << "timeStepCounter = " << timeStepCounter_ << endl;
        
        if(timeStepCounter_ >= noSnapShots_)
        {
            Info << "shifting " << endl;
            
            // shifting
            for(label i = 0; i < noSnapShots-1; ++i)
            {
                for(label j = 0; j < nMicro_+1; ++j)
                {
                    massFlow_[i][j]=massFlow_[i+1][j];
                }
            }
            
            for(label j = 0; j < nMicro_+1; ++j)
            {
                massFlow_[noSnapShots-1][j]=newM[j];
            }
            
            POD
            (
                massFlow_,
                noSnapShots_,
                nMicro_,
                newS
            );
            
            // set smoother mass flow rates due to time smoothing
            Info << "set new mass flow rates" << endl;
            
            for(label j = 0; j < nMicro_; ++j)
            {
                m_[j] = massFlowSmoother_[noSnapShots-1][j];
            }
        }
        else
        {
            for(label j = 0; j < nMicro_+1; ++j)
            {
                massFlow_[timeStepCounter_][j]=newM[j];
            }        
        }
        
        timeStepCounter_++;
    }
    
    // old part

    label M = mFourierCoeffs_.size();
    
    scalarField newS(nMicro_+1, 0.0);
    scalarField newM(nMicro_+1, 0.0); 
    
    forAll(s_, i)
    {
        newS[i] = s_[i];
        newM[i] = m_[i];
    }
    
    newS[nMicro_]=macroLength_;
    newM[nMicro_]=newM[0];        
    
    // apply fourier polynomial least squares fit through mass flow rate data points
    fourierPolyLeastSquaresFit fit(newS, newM, M, macroLength_);
    
    const scalarField& coeffs = fit.coeffs();

    forAll(coeffs, j)
    {
        mFourierCoeffs_[j] = coeffs[j];
    }
    
    Info << "mass flow rate before = " << m_ << endl;
    
    // evaluate new 'smoother' mass flow rate at each collocation node
    
    m_ = 0.0;

    forAll(s_, i)
    {
        m_[i] = fourierEvaluation(mFourierCoeffs_, i);       
    }
    
    Info << "mass flow rate after = " << m_ << endl;         
}*/


/*
scalar agentMicroElements::fourierEvaluation
(
    const scalarField& coeffs,
    const label& i
)
{
    label M = coeffs.size();
    
    label N = (M - 1)*0.5;
    
    scalar value = coeffs[0];
    
    for (label j=1; j<= N; j++)
    {
        value += coeffs[(2*j)-1]*Foam::sin(2.0*mathematicalConstant::pi*s_[i]*j/macroLength_);
        value += coeffs[(2*j)+1-1]*Foam::cos(2.0*mathematicalConstant::pi*s_[i]*j/macroLength_);                
    }
    
    return value;
}

scalar agentMicroElements::fourierDifferentiation
(
    const scalarField& coeffs,
    const label& i
)
{
    label M = coeffs.size();
    
    label N = (M - 1)*0.5;
    
    scalar diff = 0.0;
    
    for (label j=1; j<= N; j++)
    {
        diff += coeffs[(2*j)-1]*j*Foam::cos(2.0*mathematicalConstant::pi*s_[i]*j/macroLength_);
        diff -= coeffs[(2*j)+1-1]*j*Foam::sin(2.0*mathematicalConstant::pi*s_[i]*j/macroLength_);                
    }
    
    diff *= 2.0*mathematicalConstant::pi/macroLength_;
    
    return diff;
}*/

// void applyMassFlowRatePolynomialSmoothing()    
// {
// 
//         polynomialLeastSquaresFit fit(s_, m_, M);
//     
//         const scalarField& coeffs = fit.coeffs();
// 
//         forAll(coeffs, j)
//         {
//             polynomialCoeffs_[j] = coeffs[j];
//         }

//         m_ = 0.0;
//         
//         forAll(s_, i)
//         {
//             forAll(polynomialCoeffs_, j)
//             {
//                 m_[i] += (polynomialCoeffs_[j]*Foam::pow(s_[i], j));
//             }
//         }
//     }
// }
    
// periodic implementation

// void agentMicroElements::computeDensity()
// {
//     // first solve for mass flux coefficients 
//     if(fourier_) 
//     {
//         // mass flux
//         computeFourierCofficients(m_, mCoeffs_);        
//         
//         // velocity 
//         scalarField v(nMicro_, 0.0);
//         
//         forAll(v_, i)
//         {
//             v[i] = v_[i].x();
//         }
//         
//         computeFourierCofficients(v, vCoeffs_);
//     }
//     
//     if(fourierPolySmooth_) 
//     {
//         applyFourierSmoothingToMassFlux();        
//     }
//     
//     forAll(rhoM_, i)
//     {
//         if(three_)// three point stencil model
//         {
//             const scalar& m_1 = m_[iMinus(i,1)];
//             const scalar& m1 = m_[iPlus(i,1)];
//             
//             flux_[i] = -(macroDeltaT_*nReal_/(A_[i]*deltaS_))*
//                         (
//                             -m_1 + m1 
//                         )/2.0;        
//         }        
//         else if(five_) // five point stencil model
//         {
//             const scalar& m_2 = m_[iMinus(i,2)];
//             const scalar& m_1 = m_[iMinus(i,1)];
//             const scalar& m1 = m_[iPlus(i,1)];
//             const scalar& m2 = m_[iPlus(i,2)];       
//             
//             flux_[i] = -(macroDeltaT_*nReal_/(A_[i]*deltaS_))*
//                         (
//                             m_2 - 8*m_1 + 8*m1 - m2
//                         )/12.0;
//         }        
//         else if(seven_) // seven point stencil model
//         {
//             const scalar& m_3 = m_[iMinus(i,3)];        
//             const scalar& m_2 = m_[iMinus(i,2)];
//             const scalar& m_1 = m_[iMinus(i,1)];
//             const scalar& m1 = m_[iPlus(i,1)];
//             const scalar& m2 = m_[iPlus(i,2)];       
//             const scalar& m3 = m_[iPlus(i,3)];               
//             
//             flux_[i] = -(macroDeltaT_*nReal_/(A_[i]*deltaS_))*
//                         (
//                             -m_3 + 9.0*m_2 - 45*m_1 
//                             + 45*m1 - 9*m2  + m3 
//                         )/60;
//         }
//         else if(nine_) // nine point stencil model
//         {
//             const scalar& m_4 = m_[iMinus(i,4)];            
//             const scalar& m_3 = m_[iMinus(i,3)];        
//             const scalar& m_2 = m_[iMinus(i,2)];
//             const scalar& m_1 = m_[iMinus(i,1)];
//             const scalar& m1 = m_[iPlus(i,1)];
//             const scalar& m2 = m_[iPlus(i,2)];       
//             const scalar& m3 = m_[iPlus(i,3)];               
//             const scalar& m4 = m_[iPlus(i,4)];
//             
//             flux_[i] = -(macroDeltaT_*nReal_/(A_[i]*deltaS_))*
//                         (
//                             3*m_4 - 32*m_3 + 168*m_2 - 672*m_1 
//                             + 672*m1 - 168*m2  + 32*m3 - 3*m4
//                          )/840.0;
//         }
//         else if(fourier_) // fourier collocation
//         {
//             scalar diffMassFlux = 0.0;
//             
//             label N = (nMicro_-1)*0.5;
//             
//             for (label j=1; j<= N; j++)
//             {
//                 diffMassFlux += mCoeffs_[(2*j)-1]*j*Foam::cos(2.0*mathematicalConstant::pi*s_[i]*j/macroLength_);
//                 diffMassFlux -= mCoeffs_[(2*j)+1-1]*j*Foam::sin(2.0*mathematicalConstant::pi*s_[i]*j/macroLength_);                
//             }
//             
//             flux_[i] = -( macroDeltaT_*nReal_/ A_[i] )*diffMassFlux*2.0*mathematicalConstant::pi/macroLength_;
//         }
//         else if(polynomial_) // polynomial collocation (non periodicity)
//         {
// //             scalar diffMassFlux = 0.0;
// //             
// //             for (label j=1; j<nMicro_; j++)
// //             {            
// //                 diffMassFlux += scalar(j)*mCoeffs_[j]*pow(s_[i], j-1);
// //             }
// //             
// //             flux_[i] = -( macroDeltaT_*nReal_/ A_[i] )*diffMassFlux;
// 
//         }
//         
//         else if(fourierPolySmooth_)
//         {
//             scalar diffMassFlux = 0.0;
//             
//             if(euler_)
//             {
//                 diffMassFlux = fourierDifferentiation(mFourierCoeffs_, i);
//             }
// //             else if (adamsBashforth_)
// //             {
// //                 diffMassFlux = 
// //                     1.5*fourierDifferentiation(mFourierCoeffs_, i) -
// //                     0.5*fourierDifferentiation(mOldFourierCoeffs_, i);
// //             }
//             else if (adamsBashforth_)
//             {
//                 diffMassFlux = 
//                     (23.0/12.0)*fourierDifferentiation(mFourierCoeffs_, i) 
//                     - (4.0/3.0)*fourierDifferentiation(mOldFourierCoeffs_, i)
//                     + (5.0/12.0)*fourierDifferentiation(mOld2FourierCoeffs_, i);
//             }
//             
//             flux_[i] = -( diffMassFlux*macroDeltaT_*nReal_/ A_[i] );             
//         }
//         
//         else if(polynomialSmoothing_) 
//         {
//             scalar diffMassFlux = 0.0;
//             
//             for(label j=1; j <= polynomialOrder_; j++)
//             {
//                 diffMassFlux += j*polynomialCoeffs_[j]*Foam::pow(s_[i], j-1);
//             }
//             
//             flux_[i] = -( diffMassFlux*macroDeltaT_*nReal_/ A_[i] );
//         }
//         /*
//         // change in mols to insert/delete (real number)
//         
//         deltaN_[i] = (flux_[i]*bbVol_[i]/mass_) + molFluxResidual_[i];
//         
//         // NINT -> change to integer
//         if(deltaN_[i] > 0.0)
//         {
//             intDeltaN_[i] = label(deltaN_[i] + 0.5);
//         }
//         else if(deltaN_[i] < 0.0)
//         {
//             intDeltaN_[i] = label(deltaN_[i] - 0.5);        
//         }
//         
//         // account for flux residual
//         molFluxResidual_[i] = deltaN_[i] - intDeltaN_[i];
//         
// //        Info << "mol Flux = " << (flux_[i]*bbVol_[i]/mass_)
// //              << ", deltaN (with old flux residual) = " << deltaN_[i]
// //              << ", intDeltaN = " << intDeltaN_[i]
// //              << ", new mol flux residual = " << molFluxResidual_[i]
// //              << endl;
//              
//         rhoN_[i] += (intDeltaN_[i]/bbVol_[i]);
//         
//         rhoM_[i] += (mass_*intDeltaN_[i]/bbVol_[i]);
//         */
//         
//         rhoN_[i] += (flux_[i]/mass_);
//         rhoM_[i] += flux_[i];
//     }
//     
//     if(fourier_) 
//     {
//         computeFourierCofficients(rhoN_, rhoNcoeffs_);        
//     }
//     
//     if(fourierPolySmooth_)
//     {
//         computeFourierCofficients(rhoN_, rhoNFourierCoeffs_);     
//     }
// //     computeCofficients(rhoN_, rhoNcoeffs_);
//     
// //     updateHeights();
//     
// //     Info << "new densities = " << rhoN_ << endl;
// }

// void agentMicroElements::computeFourierCofficients
// (
//     const scalarField& m,
//     scalarField& mCoeffs
// )
// {
//     simpleMatrix<scalar> luMatrix(nMicro_);
//     
//     forAll(m, i)
//     {
//         luMatrix[i][0] = 1.0;            
//         
//         label N = (nMicro_-1)*0.5;
//         
//         for (label j=1; j<= N; j++)
//         {
//             luMatrix[i][(2*j) - 1] = Foam::sin(2.0*mathematicalConstant::pi*s_[i]*j/macroLength_);
//             luMatrix[i][(2*j)+1-1] = Foam::cos(2.0*mathematicalConstant::pi*s_[i]*j/macroLength_);                
//         }
//         
//         luMatrix.source()[i] = m[i];
//     }
//         
//     mCoeffs = 0.0;
// 
//     mCoeffs = luMatrix.LUsolve();          
// }

// // compute new pressures
// void agentMicroElements::computePressures()
// {
//     
// //     label M = 3;
// //     scalarField coeffs(M+1,0.0);
//     
//     // using eqn of state from channel elements
//     
//     forAll(p_, i)
//     {
//         if((rhoN_[i] >= rhoStart_[i]) && (rhoN_[i] <= rhoEnd_[i]))
//         {
//             p_[i] = 0.0;
//             
//             const scalarField& coeffs = eqnOfStateCoeffs_[i];
//             
//             label M = coeffs.size() - 1;
//             
//             forAll(coeffs, j)
//             {
//                 p_[i] += coeffs[j]*pow(rhoN_[i], (M-j));
//             }
//         }
//         else
//         {
//             FatalErrorIn("agentMicroElements::computePressures()")
//                 << " Simulation exceeded density bounds of equation of state."
//                 << nl << "Measured density = " << rhoN_[i] 
//                 << " bounds: bottom = " << rhoStart_[i]
//                 << ", top = " << rhoEnd_[i]
//                 << nl << " Stopping simualation "
//                 << exit(FatalError);         
//         }   
//     }
// }


// void agentMicroElements::computePressureGradientCorrections()
// {
//     
//     if(fourier_) 
//     {
//         computeFourierCofficients(p_, pCoeffs_);        
//     }    
//     
//     if(fourierPolySmooth_)
//     {
//         computeFourierCofficients(p_, pFourierCoeffs_);
//     }
//     
//     forAll(p_, i)
//     {
//         if(three_)
//         {
//             const scalar& p_1 = p_[iMinus(i,1)];
//             const scalar& p1 = p_[iPlus(i,1)];
//            
//             phi_[i] = -(1.0/deltaS_)*
//                         (
//                             -p_1 + p1 
//                         )/2.0;        
//         }
//         else if(five_)// five point stencil model
//         {
//             const scalar& p_2 = p_[iMinus(i,2)];
//             const scalar& p_1 = p_[iMinus(i,1)];
//             const scalar& p1 = p_[iPlus(i,1)];
//             const scalar& p2 = p_[iPlus(i,2)];       
//            
//             phi_[i] = -(1.0/deltaS_)*
//                         (
//                             p_2 - 8*p_1 + 8*p1 - p2 
//                         )/12.0;
//         }        
//         else if(seven_) // seven point stencil model
//         {
//             const scalar& p_3 = p_[iMinus(i,3)];        
//             const scalar& p_2 = p_[iMinus(i,2)];
//             const scalar& p_1 = p_[iMinus(i,1)];
//             const scalar& p1 = p_[iPlus(i,1)];
//             const scalar& p2 = p_[iPlus(i,2)];       
//             const scalar& p3 = p_[iPlus(i,3)];               
//             
//             phi_[i] = -(1.0/deltaS_)*
//                         (
//                             -p_3 + 9*p_2 -45*p_1 
//                             + 45*p1 - 9*p2  + p3 
//                         )/60.0;
//         }
//         else if(nine_) // nine point stencil model
//         {
//             const scalar& p_4 = p_[iMinus(i,4)];           
//             const scalar& p_3 = p_[iMinus(i,3)];        
//             const scalar& p_2 = p_[iMinus(i,2)];
//             const scalar& p_1 = p_[iMinus(i,1)];
//             const scalar& p1 = p_[iPlus(i,1)];
//             const scalar& p2 = p_[iPlus(i,2)];       
//             const scalar& p3 = p_[iPlus(i,3)];               
//             const scalar& p4 = p_[iPlus(i,4)];
//             
//             phi_[i] = -(1.0/deltaS_)*
//                         (
//                            3*p_4 - 32*p_3 + 168*p_2 - 672*p_1 
//                             + 672*p1 - 168*p2  + 32*p3 - 3*p4
//                         )/840.0;        
//         }
//         else if(fourier_) // fourier collocation
//         {
//             scalar diffPressure = 0.0;
//             
//             label N = (nMicro_-1)*0.5;
//             
//             for (label j=1; j<= N; j++)
//             {
//                 diffPressure += pCoeffs_[(2*j)-1]*j*Foam::cos(2.0*mathematicalConstant::pi*s_[i]*j/macroLength_);
//                 diffPressure -= pCoeffs_[(2*j)+1-1]*j*Foam::sin(2.0*mathematicalConstant::pi*s_[i]*j/macroLength_);                
//             }
//             
//             phi_[i] = -diffPressure*2.0*mathematicalConstant::pi/macroLength_;
//         }        
//         else if(polynomial_) // polynomial collocation (non periodic geometries)
//         {
// //             scalar diffPressure = 0.0;
// // 
// //             for (label j=1; j<nMicro_; j++)
// //             {            
// //                 diffPressure += scalar(j)*pCoeffs_[j]*pow(s_[i], j-1);
// //             }
// //             
// //             phi_[i] = -diffPressure;
//         }
//         else if(fourierPolySmooth_)
//         {
//             phi_[i] = -fourierDifferentiation(pFourierCoeffs_, i);
//         }
//         else if(polynomialSmoothing_) 
//         {
//             // seven point model
//             const scalar& p_3 = p_[iMinus(i,3)];        
//             const scalar& p_2 = p_[iMinus(i,2)];
//             const scalar& p_1 = p_[iMinus(i,1)];
//             const scalar& p1 = p_[iPlus(i,1)];
//             const scalar& p2 = p_[iPlus(i,2)];       
//             const scalar& p3 = p_[iPlus(i,3)];               
//             
//             phi_[i] = -(1.0/deltaS_)*
//                         (
//                             -p_3 + 9*p_2 -45*p_1 
//                             + 45*p1 - 9*p2  + p3 
//                         )/60.0;
//         }        
//     }
//     
//     if(fourier_) 
//     {
//         computeFourierCofficients(phi_, phiCoeffs_);        
//     }
//     
// //     computeCofficients(phi_, phiCoeffs_);
//         
// //     Info << "new phi = " << phi_ << endl;      
// }

// void agentMicroElements::outputPeKe(const label& i)
// {
//     if(outputEnergies_)
//     {
//         scalar kE = 0.0;
//         scalar pE = 0.0;
//         scalar mols = 0.0;    
//         
//         vector velocity = getMeanVelocity(i);
//         
//         {
//             IDLList<atomisticMolecule>::iterator mol(molecules_[i]->begin());
//         
//             for (mol = molecules_[i]->begin(); mol != molecules_[i]->end(); ++mol)
//             {
//                 if(findIndex(liquidMolIds_, mol().id()) != -1)
//                 {            
//                     const atomisticMolecule::constantProperties& constProp 
//                                     = molecules_[i]->constProps(mol().id());
//                                     
//                     mols += 1.0;
//                     pE += mol().potentialEnergy();
//                     kE += 0.5*constProp.mass()*magSqr(mol().v() - velocity);
//                 }
//             }  
//         }    
//         
//         scalar averagePe = 0.0;
//         scalar averageKe = 0.0;
//         
//         if(mols > 0.0)
//         {
//             averagePe = pE/mols;
//             averageKe = kE/mols;
//         }
// 
//         Pout << "pE = " << averagePe << ", kE = " << averageKe << endl;
//     }
// }
// 
// // get mean velocity
// vector agentMicroElements::getMeanVelocity(const label& i)
// {
//     vector mom = vector::zero;
//     scalar mass = 0.0;
//     
//     {
//         IDLList<atomisticMolecule>::iterator mol(molecules_[i]->begin());
//     
//         for (mol = molecules_[i]->begin(); mol != molecules_[i]->end(); ++mol)
//         {
//             if(findIndex(liquidMolIds_, mol().id()) != -1)
//             {
//                 const atomisticMolecule::constantProperties& constProp 
//                                 = molecules_[i]->constProps(mol().id());
//                                 
//                 mass += constProp.mass();
//                 mom += mol().v()*constProp.mass();
//             }
//         }  
//     }    
//     
//     vector meanVelocity = vector::zero;
//     
//     if(mass > 0.0)
//     {
//         meanVelocity = mom/mass;
//     }
// 
//     return meanVelocity;
// }
// 
// void agentMicroElements::controlDensity(const label& i)
// {
//     if(usher_)
//     {
//         controlDensityUsher(i);
//     }
//     else if(fade_)
//     {
//         controlDensityFade(i);
//     }
// }
// 
// void agentMicroElements::controlDensityUsher(const label& i)
// {
//     scalar mols = 0.0;
//     scalar rhoN = 0.0;
//     scalar rhoM = 0.0;
//     
//     measureDensity(mols, rhoN, rhoM, i);    
//     
//     label molsToControl = (rhoN_[i] - rhoN)*bbVol_[i];
//    
//     if(mag(molsToControl) > 0)
//     {
// /*        Pout << "mols to control = " << molsToControl << endl;   */     
//         
//         // measure velocity first (use this to set newly inserted molecules)
// //         vector meanVelocity = vector::zero;
//         
// //         if(molsToControl > 0)
// //         {
// //             meanVelocity = getMeanVelocity(i);
// //         }
//         
//         insertionScheme_[i].controlMoleculesInBoundBox
//         (
//             molecules_[i](),
//             molsToControl,
//             insertionBoxes_[i]
//         );
// 
//         // apply velocity to inserted molecules
// //         if(molsToControl > 0)        
// //         {
// //             const DynamicList<label>& tnS = insertionScheme_[i].insertedMolecules();
// //             
// //             IDLList<atomisticMolecule>::iterator mol(molecules_[i]->begin());
// //         
// //             for (mol = molecules_[i]->begin(); mol != molecules_[i]->end(); ++mol)
// //             {
// //                 if(findIndex(tnS, mol().trackingNumber()) != -1)
// //                 {
// // //                     Info << "modified molecule velocity = " <<  mol().v() << endl;
// //                     
// //                     mol().v() = meanVelocity;
// //                     
// // //                     Info << " to velocity = " <<  mol().v() << endl;
// //                 }
// //             }
// //         }
//         
//         // update molsToControl in next time-step
//         if(molsToControl > 0)
//         {
//             label nMissed = molsToControl - insertionScheme_[i].nInserted();
// //             molsMissed_[i].append(molsToControl);  
// 
// //             intDeltaN_[i] = nMissed;
//             scalar efficiency = ( scalar(molsToControl - nMissed)/scalar(molsToControl) )*100.0;
//            
//             Pout<< "[ADD] Mols to control = " << molsToControl
//                 << ", Mols inserted = " <<  insertionScheme_[i].nInserted()
//                 << ", USHER efficiency = " << efficiency
//                 << endl;
//                 
//             counterUsher_ += 1.0;
//             averageEfficiency_ += efficiency;
//             
// 
//             if(Pstream::master())
//             {
//                 scalarField timeField(1, counterUsher_);
//                 scalarField efficiencyField(1, efficiency);
//             
//                 writeTimeData
//                 (
//                     casePath_,
//                     "usher_efficiency.xy",
//                     timeField,
//                     efficiencyField,
//                     true
//                 );
//             }
//         }
//         else if(molsToControl < 0)
//         {
//             label nMissed = molsToControl + insertionScheme_[i].nDeleted();
//             
//              scalar efficiency = ( scalar(mag(molsToControl) - nMissed)/scalar(molsToControl) )*100.0;
//             
//              Pout<< "[DEL] Mols to control = " << molsToControl
//                 << ", Mols deleted = " <<  insertionScheme_[i].nDeleted()
//                 << ", USHER efficiency = " << efficiency
//                 << endl;           
//         }
//         
// //         Pout << "mols to control (next time-step) = " <<  intDeltaN_[i] << endl;          
//     }
//     
//     if(counterUsher_ >= 1.0)
//     {
//         Pout << "Average USHER efficiency = " <<    averageEfficiency_/counterUsher_ << endl;
//     }
//     
// }

// void agentMicroElements::controlDensityFade(const label& i)
// {
// 
//     fadeInsertionScheme_[i].checkFractions
//     (
//         molecules_[i]()
//     );    
//     
// //     label molsToControl = intDeltaN_[i];
// 
//     scalar mols = 0.0;
//     scalar rhoN = 0.0;
//     scalar rhoM = 0.0;
//     
//     measureDensity(mols, rhoN, rhoM, i);    
//     
//     label molsToControl = (rhoN_[i] - rhoN)*bbVol_[i];
//     
//     if(mag(molsToControl) > 0)
//     {
//         Pout << "mols to control = " << molsToControl << endl;        
//         
//         fadeInsertionScheme_[i].controlMoleculesInBoundBox
//         (
//             molecules_[i](),
//             molsToControl,
//             insertionBoxes_[i]
//         );
//         
//         // update molsToControl in next time-step
//         if(molsToControl > 0)
//         {
//             label nMissed = molsToControl - fadeInsertionScheme_[i].nInserted();
// //             intDeltaN_[i] = nMissed;
//             scalar efficiency = ( ( molsToControl - nMissed )/molsToControl )*100.0;
//            
//             Pout<< "Mols to control = " << molsToControl 
//                 << ", Mols inserted = " <<  fadeInsertionScheme_[i].nInserted()
//                 << ", FADE efficiency = " << efficiency
//                 << endl;
//                 
//             counterUsher_ += 1.0;
//             averageEfficiency_ += efficiency;                
//         }
// //         else if(molsToControl < 0)
// //         {
// //             label nMissed = molsToControl + fadeInsertionScheme_[i].nDeleted();
// //             intDeltaN_[i] = nMissed;            
// //         }
//         
// //         Pout << "mols to control (next time-step) = " <<  intDeltaN_[i] << endl;          
//     }
//     
//     if(counterUsher_ >= 1.0)
//     {
//         Pout << "Average USHER efficiency = " <<    averageEfficiency_/counterUsher_ << endl;
//     }    
// }
// 
// void agentMicroElements::controlDensities()
// {
//     if(Pstream::parRun())
//     {
//         label p = Pstream::myProcNo();
//         
//         forAll(procsMDaddressing_[p], j)
//         {
//             label i = procsMDaddressing_[p][j];
//             
//             controlDensity(i);
//         }    
//     }
//     else
//     {
//         forAll(molecules_, i)
//         {
//             controlDensity(i);
//         }
//     }
// }


// FORCES

// void agentMicroElements::computeExternalForce(const label& i)
// {
//     gravityModel_->updateForce();
// }
// 
// void agentMicroElements::computeExternalForces()
// {
//     if(Pstream::parRun())
//     {
//         label p = Pstream::myProcNo();
//         
//         forAll(procsMDaddressing_[p], j)
//         {
//             label i = procsMDaddressing_[p][j];
//             
//             computeExternalForce(i);
//         }    
//     }
//     else
//     {
//         forAll(molecules_, i)
//         {
//             computeExternalForce(i);
//         }
//     } 
// }

// void agentMicroElements::controlExternalForce(const label& i)
// {
//     // leap frog implementation - apply force only after the first half time-step
//     bool control = true;
//     
//     if(startFromZero_)
//     {
//         if(microTimes_[i] <= 0.5*N_*microDeltaT_)
//         {
//             control = false;
//         }
//     }
//     
//     if(control)
//     {
//         vector gravityForce = gravityModel_->force(macroTime_);
//         
//         Info << "Force = " << gravityForce << endl;
// 
//         f_[i] = gravityForce + ((phi_[i]/rhoN_[i]) )*unitVector_; 
// 
//         IDLList<atomisticMolecule>::iterator mol(molecules_[i]->begin());
//     
//         for (mol = molecules_[i]->begin(); mol != molecules_[i]->end(); ++mol)
//         {
//             if(findIndex(liquidMolIds_, mol().id()) != -1)
//             {
//                 const atomisticMolecule::constantProperties& constProp 
//                                 = molecules_[i]->constProps(mol().id());
//                 
//                 mol().f() += f_[i];
//                 mol().a() += f_[i]/constProp.mass();
//             }
//         }
//     }
// }

// void agentMicroElements::controlExternalForces()
// {
//     forAll(molecules_, i)
//     {
//         controlExternalForce(i);
//     }
// }
// 
// void agentMicroElements::measureMassFlux(const label& i)
// {
//     scalar mass = 0.0;
//     vector mom = vector::zero;
//     
//     IDLList<atomisticMolecule>::iterator mol(molecules_[i]->begin());
// 
//     for (mol = molecules_[i]->begin(); mol != molecules_[i]->end(); ++mol)
//     {
//         if(findIndex(liquidMolIds_, mol().id()) != -1)
//         {
//             const atomisticMolecule::constantProperties& constProp 
//                             = molecules_[i]->constProps(mol().id());
//                             
//             mass += constProp.mass();                
//             mom += mol().v()*constProp.mass();
//         }
//     }
//     
//     scalar massFlux = (mom & unitVector_)/microLength_;
//     
//     vector velocity = vector::zero;
//     
//     if(mass > 0.0)
//     {
//         velocity = mom/mass;        
//     }
//     
//     m_[i] += massFlux/((N_+1)*nReal_);
//     v_[i] += velocity/(N_+1);
//     
// //     Info << "mass flux = " << m_[i] << ", vel = " << v_[i] << endl;
// }
/*
void agentMicroElements::measureMassFluxes()
{
    if(Pstream::parRun())
    {
        label p = Pstream::myProcNo();
        
        forAll(procsMDaddressing_[p], j)
        {
            label i = procsMDaddressing_[p][j];
            
            measureMassFlux(i);
        }    
    }
    else
    {
        forAll(molecules_, i)
        {
            measureMassFlux(i);
        }
    }
}*/

// void agentMicroElements::pseudoMassFlux()
// {
//     forAll(m_, i)
//     {
//         scalar random = molecules_[i]->rndGen().GaussNormal();
//         random = 0.0;
//         mOld_[i] = mOld_[i] + microDeltaT_*((-mOld_[i]/A_[i]) + mag(f_[i])) + random;
//         
//         m_[i] += mOld_[i]/(N_+1);        
//         
//         Info << " massflux (instant) = " << mOld_[i]  << ", average = " << m_[i] << endl;
//     }
// }

// void agentMicroElements::computeDensityProperty(const label& i)
// {
//     scalar rAv = 0.0;
//     scalar rMin = GREAT;
//     scalar mols = 0.0;
//     scalar pE = 0.0;
//     scalar mass = 0.0;
//     
//     vector momentum = vector::zero;
//     
//     IDLList<atomisticMolecule>::iterator mol(molecules_[i]->begin());
// 
//     for (mol = molecules_[i]->begin(); mol != molecules_[i]->end(); ++mol)
//     {
//         if(insertionBoxes_[i].contains(mol().position()))
//         {
//             const atomisticMolecule::constantProperties& constProp 
//                                 = molecules_[i]->constProps(molId_);
//                                 
//             if(mol().fraction() == 1.0)
//             { 
//                 mols += 1.0;                
//                 mass += constProp.mass();
//             	pE += mol().potentialEnergy();
//                 momentum += constProp.mass()*mol().v();
//                 
//                 if(mol().R() < rMin)
//                 {
//                     rMin = mol().R();
//             	}
//             
//                 rAv += mol().R();
//             }
//         }
//     }
//     
//     vector velocity = vector::zero;
//     
//     if(mols > 0.0)
//     {
//         pE /= mols;
//         rAv /= mols;
//         velocity = momentum/mass;
//     }
//     
//     
//     // USHER
//     
//     insertionScheme_[i].rChar() = rMin;
//     insertionScheme_[i].targetPE() = pE;
//     
//     scalar rho = mass/bbVol_[i];
//     insertionScheme_[i].deltaS() = 0.1/Foam::pow(rho,1.5);
//     
//     insertionScheme_[i].velocity() = velocity;     
// /*        
//     Info << "USHER information. rChar = " << insertionScheme_[i].rChar() 
//             << ", target pE " << insertionScheme_[i].targetPE()
//             << ", delta S = " << insertionScheme_[i].deltaS()
//             << endl;*/
// 
//     // FADE
//     
//     fadeInsertionScheme_[i].Rref() = rAv;
//     fadeInsertionScheme_[i].targetPE() = pE;    
// }               

// void agentMicroElements::computeDensityProperties()
// {
//     if(Pstream::parRun())
//     {
//         label p = Pstream::myProcNo();
//         
//         forAll(procsMDaddressing_[p], j)
//         {
//             label i = procsMDaddressing_[p][j];
//             
//             computeDensityProperty(i);
//         }
//     }
//     else
//     {
//         forAll(molecules_, i)
//         {
//             computeDensityProperty(i);
//         }    
//     }
// }
// 
// void agentMicroElements::evolveBeforeForce(const label& i)
// {
//     molecules_[i]->evolveBeforeForces();
// }
// 
// void agentMicroElements::evolveBeforeForces()
// {
//     forAll(molecules_, i)
//     {
//         molecules_[i]->evolveBeforeForces();
//     }
// }
// 
// 
// void agentMicroElements::updateAcceleration(const label& i)
// {
//     molecules_[i]->updateAcceleration();
// }
// 
// void agentMicroElements::updateAccelerations()
// {
//     forAll(molecules_, i)
//     {
//         molecules_[i]-> updateAcceleration();
//     }
// }
// 
// void agentMicroElements::evolveAfterForce(const label& i)
// {
//     molecules_[i]-> controlAfterForces();
//     molecules_[i]->finalHalfVelocity();
//     
//     thermostat(i);
//     
//     molecules_[i]->postPreliminaries();
//     
//     outputPeKe(i);
// }
// 
// 
// void agentMicroElements::evolveAfterForces()
// {
//     forAll(molecules_, i)
//     {
//         evolveAfterForce(i);
//     }
// }
// 
// 
// 
// void agentMicroElements::calculateForce(const label& i)
// {
//     molecules_[i]->calculateForce();  
// }
// 
// void agentMicroElements::calculateForces()
// {
//     forAll(molecules_, i)
//     {
//         molecules_[i]->calculateForce();
//     }
// }
// 
// 
// 
// void agentMicroElements::syncMassFluxes()
// {
//     
//     if(Pstream::parRun())
//     {
//         // parallel communications
//         forAll(m_, i)
//         {
//             reduce(m_[i], sumOp<scalar>());
//             reduce(v_[i], sumOp<vector>());
//             reduce(f_[i], sumOp<vector>());            
//         }    
//     }
// }
// 
// 


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const label& agentMicroElements::nMicro() const
{
    return nMicro_;
}

// const label& agentMicroElements::nRealisations() const
// {
//     return nReal_;
// }
/*
void agentMicroElements::output()
{
//     if(runTime_.outputTime())
    {
        // measure densities
        scalarField measMols(nMicro_, 0.0);
        scalarField measRhoN(nMicro_, 0.0);
        scalarField measRhoM(nMicro_, 0.0);
        scalarField measPressure(nMicro_, 0.0);
        scalarField measTemperature(nMicro_, 0.0);        
        
        measureDensities(measMols, measRhoN, measRhoM);
        thermoMeasurements(measTemperature, measPressure);
       
        
        if(Pstream::master())
        {
            forAll(molecules_, i)
            {
                scalarField massFlux(1, m_[i]);
                scalarField pressure(1, p_[i]);
                scalarField rhoN(1, rhoN_[i]);
                scalarField rhoM(1, rhoM_[i]);
                scalarField phi(1, phi_[i]);
                scalarField f(1, f_[i].x()); 
                scalarField v(1, v_[i].x());
                
                scalarField measMolsOutput(1, measMols[i]);
                scalarField measRhoNOutput(1, measRhoN[i]);
                scalarField measRhoMOutput(1, measRhoM[i]);  
                scalarField measPressureOutput(1, measPressure[i]); 
                scalarField measTempOutput(1, measTemperature[i]);                 

                
                scalarField timeField(1, macroTime_);
                
                std::string s;
                std::stringstream out;
                out << i;
                s = out.str(); 
                
                writeTimeData
                (
                    casePath_,
                    "massFlux_"+s+".xy",
                    timeField,
                    massFlux,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "pressure_"+s+".xy",
                    timeField,
                    pressure,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "rhoN_"+s+".xy",
                    timeField,
                    rhoN,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "rhoM_"+s+".xy",
                    timeField,
                    rhoM,
                    true
                );       
                
                writeTimeData
                (
                    casePath_,
                    "phi_"+s+".xy",
                    timeField,
                    phi,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "force_"+s+".xy",
                    timeField,
                    f,
                    true
                );  
                
                writeTimeData
                (
                    casePath_,
                    "velocity_"+s+".xy",
                    timeField,
                    v,
                    true
                );
                

                // measured properties
                
                writeTimeData
                (
                    casePath_,
                    "meas_mols_"+s+".xy",
                    timeField,
                    measMolsOutput,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "meas_rhoN_"+s+".xy",
                    timeField,
                    measRhoNOutput,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "meas_rhoM_"+s+".xy",
                    timeField,
                    measRhoMOutput,
                    true
                );                
                writeTimeData
                (
                    casePath_,
                    "meas_p_"+s+".xy",
                    timeField,
                    measPressureOutput,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "meas_T_"+s+".xy",
                    timeField,
                    measTempOutput,
                    true
                );                
            }
            
            if(fourierPolySmooth_)
            {
                // output coefficients
                forAll(mFourierCoeffs_, i)
                {  
                    scalarField mCoeffs(1, mFourierCoeffs_[i]);
                    
                    scalarField timeField(1, macroTime_);
                    
                    std::string s;
                    std::stringstream out;
                    out << i;
                    s = out.str();                 
                    
                    
                    writeTimeData
                    (
                        casePath_,
                        "massFlux_coeffs_"+s+".xy",
                        timeField,
                        mCoeffs,
                        true
                    );                  
                }
                
                forAll(pFourierCoeffs_, i)
                {  
                    scalarField pCoeffs(1, pFourierCoeffs_[i]);
                    scalarField rhoNCoeffs(1, rhoNFourierCoeffs_[i]);                    
                    
                    scalarField timeField(1, macroTime_);
                    
                    std::string s;
                    std::stringstream out;
                    out << i;
                    s = out.str();                 
                    
                    
                    writeTimeData
                    (
                        casePath_,
                        "pressure_coeffs_"+s+".xy",
                        timeField,
                        pCoeffs,
                        true
                    );                  
                    
                    writeTimeData
                    (
                        casePath_,
                        "rhoN_coeffs_"+s+".xy",
                        timeField,
                        rhoNCoeffs,
                        true
                    );                     
                }                
               
            }
            
            // coefficients 
            if(fourier_)
            {
                forAll(mCoeffs_, i)
                {            
                    scalarField mCoeffs(1, mCoeffs_[i]);
                    scalarField pCoeffs(1, pCoeffs_[i]);                
                    scalarField vCoeffs(1, vCoeffs_[i]);                 
                    scalarField rhoNCoeffs(1, rhoNcoeffs_[i]);                 
                    scalarField phiCoeffs(1, phiCoeffs_[i]); 

                    scalarField timeField(1, macroTime_);
                    
                    std::string s;
                    std::stringstream out;
                    out << i;
                    s = out.str();                 
                    
                    
                    writeTimeData
                    (
                        casePath_,
                        "massFlux_coeffs_"+s+".xy",
                        timeField,
                        mCoeffs,
                        true
                    );
                    
                    writeTimeData
                    (
                        casePath_,
                        "velocity_coeffs_"+s+".xy",
                        timeField,
                        vCoeffs,
                        true
                    );
                    
                    writeTimeData
                    (
                        casePath_,
                        "rhoN_coeffs_"+s+".xy",
                        timeField,
                        rhoNCoeffs,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        "pressure_coeffs_"+s+".xy",
                        timeField,
                        pCoeffs,
                        true
                    );                

                    writeTimeData
                    (
                        casePath_,
                        "phi_coeffs_"+s+".xy",
                        timeField,
                        phiCoeffs,
                        true
                    );                 
                }
            }
            
            scalarField M(1, 0.0); 
            
            forAll(rhoM_, i)
            {
                M[0] += A_[i]*rhoM_[i];
            }
            
            M[0] *= macroLength_/nMicro_;
            
            scalarField timeField(1, macroTime_);  
                
            writeTimeData
            (
                casePath_,
                "M.xy",
                timeField,
                M,
                true
            );     
            
            
            // CAI OUTPUTTING
            
            if(CAI_)
            {
                scalarField S(1, S_);
                scalarField Tmacro(1, Tmacro_);
                scalarField period(1, Tmacro_);
                scalarField gMicro(1, gMicro_);
//                 scalarField force(1, mag(caiForcing_));
                scalarField N(1, N_);
                scalarField macroDeltaT(1, macroDeltaT_);
                
                scalarField timeField(1, macroTime_);

                writeTimeData
                (
                    casePath_,
                    "scaleSeparation.xy",
                    timeField,
                    S,
                    true
                ); 
                
                writeTimeData
                (
                    casePath_,
                    "Tmacro.xy",
                    timeField,
                    Tmacro,
                    true
                ); 
                
                writeTimeData
                (
                    casePath_,
                    "period.xy",
                    timeField,
                    period,
                    true
                ); 
                
                writeTimeData
                (
                    casePath_,
                    "gMicro.xy",
                    timeField,
                    gMicro,
                    true
                ); 
                
//                 writeTimeData
//                 (
//                     casePath_,
//                     "force.xy",
//                     timeField,
//                     force,
//                     true
//                 ); 
                
                writeTimeData
                (
                    casePath_,
                    "N.xy",
                    timeField,
                    N,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "macroDeltaT.xy",
                    timeField,
                    macroDeltaT,
                    true
                );            
            }
        }
    }
}

*/

/*
void agentMicroElements::clear()
{
//     mOld_ = m_;

    mOld2FourierCoeffs_ = mOldFourierCoeffs_;
    mOldFourierCoeffs_ = mFourierCoeffs_;
    
    m_ = 0.0;
    v_ = vector::zero;
    f_ = vector::zero;    
//     T_ = 0.0;
}*/




// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
// measure density of all micro elements
void agentMicroElements::measureDensities
(
    scalarField& mols,
    scalarField& rhoN,
    scalarField& rhoM
)
{
    if(Pstream::parRun())
    {
        mols = 0.0;
        rhoN = 0.0;
        rhoM = 0.0;
        
        label p = Pstream::myProcNo();
        
        forAll(procsMDaddressing_[p], j)
        {
            label i = procsMDaddressing_[p][j];
            
            measureDensity
            (
                mols[i],
                rhoN[i],
                rhoM[i],
                i
            );
        }
        
        // parallel communications
        forAll(mols, i)
        {
            reduce(mols[i], sumOp<scalar>());
            reduce(rhoN[i], sumOp<scalar>());
            reduce(rhoM[i], sumOp<scalar>());
        }
    }
    else
    {
        forAll(molecules_, i)
        {
            measureDensity
            (
                mols[i],
                rhoN[i],
                rhoM[i],
                i
            );
        }
    }
    
//    Info<< "Measuring density: nMols = " << mols
//        << " rhoN = " << rhoN
 //       << " rhoM = " << rhoM        
  //      << endl; 
}*/

// measure density of one micro element
/*

void agentMicroElements::measureDensity
(
    scalar& mols,
    scalar& rhoN,
    scalar& rhoM,
    const label& i
)
{
    label nMols = 0;
    scalar mass = 0.0;       
    
    IDLList<atomisticMolecule>::iterator mol(molecules_[i]->begin());

    for (mol = molecules_[i]->begin(); mol != molecules_[i]->end(); ++mol)
    {
        if(findIndex(liquidMolIds_, mol().id()) != -1)
        {
            const atomisticMolecule::constantProperties& constProp 
                            = molecules_[i]->constProps(mol().id());
            nMols++;
            
            mass += constProp.mass();
        }
    }
    
    mols = nMols;
    rhoN = nMols/bbVol_[i];
    rhoM = mass/bbVol_[i];
}*/

/*
void agentMicroElements::thermoMeasurements
(
    scalarField& temperatureField,
    scalarField& pressureField
)
{
    if(Pstream::parRun())
    {
        temperatureField = 0.0;
        pressureField = 0.0;        
        
        label p = Pstream::myProcNo();
        
        forAll(procsMDaddressing_[p], j)
        {
            label i = procsMDaddressing_[p][j];
            
            measurePressAndTemp
            (
                temperatureField[i],            
                pressureField[i],
                i
            );
        }
        
        // parallel communications
        forAll(pressureField, i)
        {
            reduce(temperatureField[i], sumOp<scalar>());
            reduce(pressureField[i], sumOp<scalar>());            
        }
    }
    else
    {
        forAll(molecules_, i)
        {
            measurePressAndTemp
            (
                temperatureField[i],            
                pressureField[i],
                i
            );
        }
    }
    
    Info<< "Measuring Pressure: " << pressureField
        << endl;  
        
    Info<< "Measuring temperature: " << temperatureField
        << endl;        
}*/

/*
void agentMicroElements::measurePressAndTemp
(
    scalar& temperature,
    scalar& pressure,
    const label& i
)
{
    vector velocity = vector::zero;
    
    {
        scalar mass = 0.0;
        vector mom = vector::zero;
        
        IDLList<atomisticMolecule>::iterator mol(molecules_[i]->begin());

        for (mol = molecules_[i]->begin(); mol != molecules_[i]->end(); ++mol)
        {
            if(findIndex(liquidMolIds_, mol().id()) != -1)
            {
                const atomisticMolecule::constantProperties& constProp 
                                = molecules_[i]->constProps(mol().id());
                                
                mass += constProp.mass();
                mom += constProp.mass()*mol().v();
            }
        }
        
        velocity = mom/mass;
    }
    
    label nMols = 0;
    scalar mass = 0.0;  
    scalar dof = 0.0;
    scalar kE = 0.0;    
    tensor kineticTensor = tensor::zero;
    tensor virialTensor = tensor::zero;
    
    IDLList<atomisticMolecule>::iterator mol(molecules_[i]->begin());

    for (mol = molecules_[i]->begin(); mol != molecules_[i]->end(); ++mol)
    {
        if(findIndex(liquidMolIds_, mol().id()) != -1)
        {
            const atomisticMolecule::constantProperties& constProp 
                            = molecules_[i]->constProps(mol().id());
            nMols++;
            
            mass += constProp.mass();
            
            dof += constProp.degreesOfFreedom();
            
            kE += 0.5*constProp.mass()*magSqr(mol().v() - velocity);
            
            kineticTensor += ( constProp.mass()*(mol().v()-velocity)*(mol().v()-velocity) );
            virialTensor += 0.5*mol().rf();
        }
    }

    const scalar& kB = molecules_[i]->redUnits().kB();
    
    if(nMols > 0.0)
    {
        pressure = tr( (3.0*nMols*kineticTensor/dof) + virialTensor)
                                    /( 3.0*bbVol_[i]);    
                                    
        temperature = (2.0*kE)/(kB*dof);    
    }
}

*/


/*
void agentMicroElements::readInThermostatProperties()
{
    temperature_ = readScalar(propsDict_.lookup("temperature"));    
    tauT_ = readScalar(propsDict_.lookup("tauT"));
    scalar binWidth = readScalar(propsDict_.lookup("temperatureBinWidth"));
    
    
//     List<label> bins = List<label>(propsDict_.lookup("temperatureBins"));     
    
    tempBins_.setSize(nMicro_);
    tempBinWidths_.setSize(nMicro_);
    tempStartPoints_.setSize(nMicro_);    
    
//     if( bins.size() != nMicro_ )
//     {
//         FatalErrorIn("agentMicroElements::agentMicroElements()")
//             << "number of bins = " << bins.size() 
//             <<" not equal to no of micro elements = " << nMicro_
//             << exit(FatalError);
//     }

    forAll(tempBins_, i)
    {
        
        tempBins_[i] = label( (h_[i]/binWidth) + 0.5);
        
        if( tempBins_[i] == 0 )
        {
            FatalErrorIn("agentMicroElements::agentMicroElements()")
                << "number of bins = " << tempBins_[i] 
                <<" has to be greater or equal to 1. Try increasing binWidth " << binWidth
                << exit(FatalError);
        }       
        
        tempBinWidths_[i]=h_[i]/tempBins_[i];
        
        const boundBox& bb = boxes_[i];
        
        vector startPoint = bb.midpoint()-0.5*h_[i]*normal_;
        
        tempStartPoints_[i]= startPoint;
        
        Info << "start point = " << startPoint << endl;
    }
    
}*/

/*

label agentMicroElements::isPointWithinBin
(
    const vector& rI,
    const label& i
)
{
    label binNumber = -1;

    vector rSI = rI - tempStartPoints_[i];
    scalar rD = rSI & normal_;
    label n = label(rD/tempBinWidths_[i]);

    if
    (
        (n >= 0) && (rD <= h_[i])
    )
    {

        if(n == tempBins_[i]) 
        {
            n--;
        }

        if(n < tempBins_[i])
        {
            binNumber = n;
        }
    }

    return binNumber;
}*/

/*
void agentMicroElements::readInEquationOfStates()
{
    label degree = readLabel(propsDict_.lookup("eqnOfStateDegree"));
    
    PtrList<entry> infoList(propsDict_.lookup("eqnOfStates"));

    if( infoList.size() != nMicro_ )
    {
        FatalErrorIn("agentMicroElements::agentMicroElements()")
            << "Number of equations of state, = " << infoList.size()  
            <<" not equal to no of micro elements = " << nMicro_
            << exit(FatalError);
    }    
    rhoStart_.setSize(nMicro_);
    rhoEnd_.setSize(nMicro_);
    
    forAll(infoList, i)
    {
        const entry& stateI = infoList[i];
        const dictionary& dict = stateI.dict();
        
        eqnOfStateCoeffs_[i].setSize(degree, 0.0);
        
        List<scalar> coeffs = List<scalar>(dict.lookup("coeffs"));   
        
        if( coeffs.size() != degree )
        {
            FatalErrorIn("agentMicroElements::agentMicroElements()")
                << "Equation of state coeffs = " << coeffs 
                <<" not equal to degree = " << degree
                << exit(FatalError);
        }

        forAll(coeffs, j)
        {
            eqnOfStateCoeffs_[i][j] = coeffs[j];
        }
        if (dict.found("rhoStart"))
        {
            rhoStart_[i] = readScalar(dict.lookup("rhoStart"));    
            rhoEnd_[i] = readScalar(dict.lookup("rhoEnd"));         
        }
        else // default
        {
            rhoStart_[i] = 0.3;
            rhoEnd_[i] = 0.95;
        }
    }
}*/

// void agentMicroElements::setBoundBoxes()
// {
//     PtrList<entry> boxList(propsDict_.lookup("boxes"));
// 
//     boxes_.setSize(boxList.size());
//     insertionBoxes_.setSize(boxList.size());
//     
//     if( boxList.size() != nMicro_ )
//     {
//         FatalErrorIn("agentMicroElements::agentMicroElements()")
//             << "Number of boxes, = " << boxList.size()  
//             <<" not equal to no of micro elements = " << nMicro_
//             << exit(FatalError);
//     }    
//     
//     forAll(boxList, b)
//     {
//         const entry& boxI = boxList[b];
//         const dictionary& dict = boxI.dict();
// 
//         vector startPoint = dict.lookup("startPoint");
//         vector endPoint = dict.lookup("endPoint");
//         checkBoundBox(boxes_[b], startPoint, endPoint);
//     }
//     
//     
//     // set areas and volumes
//     forAll(boxes_, i)
//     { 
//         const boundBox& bb = boxes_[i];
//         
//         h_[i] = bb.span().y();
// 
//         A_[i] = h_[i] * bb.span().z(); 
//         
//         bbVol_[i] = A_[i]*microLength_;
//         
//         insertionBoxes_[i].min() = bb.min() + offset_*bb.span().y()*normal_;
//         insertionBoxes_[i].max() = bb.max() - offset_*bb.span().y()*normal_;
//     }
//         
//     
// }






void agentMicroElements::parallelProcessingRequirements()
{
    if(nMicro_ < (Pstream::nProcs() - 1))
    {
        FatalErrorIn("agentMicroElements::agentMicroElements()")
            << "There are too many processors (" << Pstream::nProcs() 
            << ") available for the number of MD subdomains (" <<  nMicro_
            << ")." << nl
            << exit(FatalError);
    }

    //round 1: compute and assign the "floor" number of MD subdomains per processor

    label nM = nMicro_;
    label nP = Pstream::nProcs();
    scalar nMperProcScalar = scalar(nM)/scalar(nP);

    label nMperProc = label(nMperProcScalar+0.5);
    
    Info << "no of micro elements per processor = " << nMperProc << endl;

    label mCounter = 0;

    // for all processors
    for (label p=0; p < nP; p++)
    {
        for (label j=0; j < nMperProc; j++)
        {
            processorNumbers_[mCounter] = p;

            mCounter++;
        }
    }

    //round 2: assign the residual number of MD subdomains to the rest of the processors

    label nMresidual = nM - (nMperProc*nP);

    for (label p=0; p < nMresidual; p++)
    {
        processorNumbers_[mCounter] = p;

        mCounter++;
    }

    Info << nl << " processorNumbers: " << processorNumbers_ << endl;


    //- building the processor-to-MDsubdomain id addressing

    List<DynamicList <label> > procsMDaddressing(nP);

    forAll(processorNumbers_, m)
    {
        procsMDaddressing[processorNumbers_[m]].append(m);
    }

    forAll(procsMDaddressing, p)
    {
        procsMDaddressing_[p].transfer(procsMDaddressing[p].shrink());
    }

    Info << "procsMDaddressing: " << procsMDaddressing_ << endl;
}
// const vectorField& agentMicroElements::positions() const
// {
//     return positions_;
// }
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
