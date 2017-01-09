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
    icoMacroModel

Description

    this is the incompressible macroscopic model used for use in the IMM

    enter all values in SI units
    

\*---------------------------------------------------------------------------*/

#ifndef icoMacroModel_H
#define icoMacroModel_H


#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "boundBox.H" 
#include "reducedUnits.H"
#include "scalarMatrices.H"
#include "SVD.H"
#include "SortableList.H"
#include "simpleMatrix.H"
#include "mathematicalConstants.H"
#include "writeTimeData.H"

#include "Random.H"
#include "polynomialLeastSquaresFit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class icoMacroModel Declaration
\*---------------------------------------------------------------------------*/

class icoMacroModel
{

    // Private data
        Time& time_;
        
        const reducedUnits& rU_;
        
	    //- subDictionary containing the properties
//         dictionary propsDict_;

        fileName casePath_;
        
        // no of micro elements
        label nMicro_;

        // number of micro element steps
        label N_;
        
        label nIter_;
        label iter_;
        
        // length of macro channel
        scalar L_;
        
        // uniform spacing between micro elements
        scalar deltaS_;

        scalar Fs_;
        
        // reference number density
        scalar rhoN_;

        // unit vector pointing in direction of flow
        vector n_;
        
        // collocation points
        scalarField s_;

        // mass flow rates
        scalarField mDot_;

        // pressure corrector terms
        scalarField phi_;

        // micro forces
        scalarField f_;
        List<vector> forces_;
        
        // coefficients
        scalarField k_;   
        
        
        
        //- OUTPUT
        
        // outerList = nIters
        scalarField mDotMean_;
        scalarField iteration_;

        List<scalarField> mDotAll_;        
        List<scalarField> phisAll_;
        List<scalarField> phisCoeffsAll_;        
        List<scalarField> fAll_;        
        
        label nPts_;
        
        scalarField sSmooth_;        
        List<scalarField> phisSmooth_;
        
    // Private Member Functions
                  

public:

    // Constructors

        //- Constructor 
        icoMacroModel
        (
            Time& t,
            const reducedUnits& rU,
            const dictionary& dict
        );
        
        ~icoMacroModel();
        
    // Member Functions
        
        
        
        void solve();
        
//         void pseudoMD();
        

        const label& nMicro() const;  
        
        const label& nIter() const;

        // to receive from MD
        scalarField& massFlowRates();
        
        // to send to MD
        List<vector>& setAndGetForces();
        
        
        // write out data
        
        void write();
        
        void writeVariables
        (
            const fileName& pathName,
            const word& nameFile,    
            const scalarField& xData,
            const List<scalarField>& yData
        );
        
        void writeVariablesSideways
        (
            const fileName& pathName,
            const word& nameFile,    
            const scalar& xData,
            const scalarField& yData
        );    
        
        

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
