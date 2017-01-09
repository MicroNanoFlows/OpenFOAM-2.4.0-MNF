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
    agentMicroElements

Description

        use for IMM
	


    

\*---------------------------------------------------------------------------*/

#ifndef agentMicroElements_H
#define agentMicroElements_H


#include "List.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "polyMesh.H"
#include "agentCloud.H"
// #include "boundBox.H" 
// #include "atomisticUsher.H"
// #include "gravityForce.H" 
// #include "atomisticFade.H"

// #include "scalarMatrices.H"
// #include "SVD.H"
// #include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class agentMicroElements Declaration
\*---------------------------------------------------------------------------*/

class agentMicroElements
{

    // Private data
        
        // MD time scheme
        Time& runTime_;

        const agentProperties& cP_;
        
        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary immDict_;        
        
        //- subDictionary containing the properties
        dictionary propsDict_;

        fileName casePath_;
        
        // no of micro elements
        label nMicro_;
        
        // List of clouds
        List< autoPtr<agentCloud> > clouds_;
        
        // number of micro element steps
        label N_;
        
        label myProcNo_;
        
        // length of macro channel
        scalar macroLength_;
        
        // macro time-step
        scalar macroDeltaT_;
        
        // when the macro solution should end
        scalar macroEndTime_;
        
        // update of time
        scalar macroTime_;
        
        scalarField microTimes_;        
        
        // reference density
//         scalar densityRef_;

        
        // micro time-step
        scalar microDeltaT_;
        
        // uniform spacing between micro elements
        scalar deltaS_;
        
        // micro element locations
        scalarField s_;
        
        // number of agents
        scalarField nAgents_;
        
        // mass density
        scalarField rhoM_; 
        
        
        // measured mass fluxes
        scalarField m_;
//         scalarField mCoeffs_;
        

//         scalarField deltaN_;
//         labelField intDeltaN_; 
//         scalarField molFluxResidual_;
        
        // cross sectional width
        scalarField W_;
        
        // liquid molecule ids
//         List<label> agentIds_;
//         label molId_;
        
        // mass of molecule(s)
//         scalar mass_;
        
        // liquid volume 
//         List<boundBox> boxes_;
//         List<scalar> bbVol_;
        
        // equation of state coefficients
//         List<scalarField> eqnOfStateCoeffs_;
        
        // micro 
        vector unitVector_;
        scalar microLength_;

        
        // list of clouds - each entry stores its processor number
        labelList processorNumbers_;

        // specify proc number, receive list of cloud ids
        labelListList procsMDaddressing_;


//         List<scalar> rhoStart_;
//         List<scalar> rhoEnd_; 
        
        
//         bool CAI_;
//         label nMacro_;
//         scalar kG_;
//         
//         // input
//         scalarField TMacro_;
//         scalarField timeMacro_;        
//         
//         scalar Tmacro_;
//         scalar Tmicro_;
//         scalar S_;
//         scalar gMicro_;
        
        
//         bool startFromZero_;
        
        
    // Private Member Functions
            
//             void setBoundBoxes();

//             void checkBoundBox
//             (
//                 boundBox& b,
//                 const vector& startPoint,
//                 const vector& endPoint
//             );

                    
            void parallelProcessingRequirements();
            
//             void setNumericalDifferentiationMethod();
            
//             void setNumericalIntegrationMethod();

//             void readInEquationOfStates();
            
//             void readInThermostatProperties();
            
//             void readVariableTimeSteppingInfo();
/*            
            label isPointWithinBin
            (
                const vector& rI,
                const label& i
            );*/
            
//             void applyFourierSmoothingToMassFlux();
            
//             void getDeltaH(const scalar& rhoN, scalar& deltaH);
/*
            void computeCofficients
            (
                const scalarField& m,
                scalarField& mCoeffs
            );*/

//             void computeFourierCofficients
//             (
//                 const scalarField& m,
//                 scalarField& mCoeffs
//             );          
//             
//             vector getMeanVelocity(const label& i);
//             
//             void outputPeKe(const label& i);            
//             
//             void setInsertionMethod();
//             
//             void controlDensityUsher(const label& i);
//             void controlDensityFade(const label& i);
// 
//             
//             scalar fourierEvaluation
//             (
//                 const scalarField& coeffs,
//                 const label& i
//             );
//             
//             scalar fourierDifferentiation
//             (
//                 const scalarField& coeffs,
//                 const label& i
//             );
//             
//             void setTMacro(); 
            


public:

    // Constructors

        //- Null Constructor 
//         agentMicroElements();

        //- Constructor 
        agentMicroElements
        (
            Time& t,
            const polyMesh& mesh,
            const agentProperties& cP         
//             const dictionary& dict
            
        );
        
        ~agentMicroElements();
        
    // Member Functions
            
//             void clear();
            
//             void output();
                        
            void updateMacroTime();
            void updateMicroTime(const label& i);
            void updateMicroTimes();
            
            
        // Access

            const label& nMicro() const;
//             const label& nRealisations() const;

//             inline const bool& pseudo() const;
            
            inline bool run ();
            
            inline const scalar& macroTime() const;
            inline const scalar& microTime(const label& i) const;
            
            inline const label& N() const;            
            
            inline label iMinus (const label& i, const label& value);
            inline label iPlus (const label& i, const label& value);

            inline List< autoPtr<agentCloud> >& clouds();
            inline const List< autoPtr<agentCloud> >& clouds() const;
            
            
            inline const labelListList& procsMDaddressing() const;    

  
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "agentMicroElementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
