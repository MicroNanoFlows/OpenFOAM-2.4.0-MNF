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
    nameFile5_(propsDict_.lookup("nameFile_referenceMols"))
{
//     endTime_ = time_.endTime().value();
    
//     Info << "endTime = " << endTime_ << endl;
    
//     deltaT_ = time_.deltaT().value();
    
//     nSteps_ = endTime_/deltaT_;
    
//     nOutputSteps_=1000;
//     
//     if(propsDict_.found("nOutputSteps"))
//     {
//         nOutputSteps_ = readLabel(propsDict_.lookup("nOutputSteps"));
//     }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dropletInfo::~dropletInfo()
{}





// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dropletInfo::createField()
{
    readFromStorage();
 
}


void dropletInfo::meanFreePathVsTime()
{
    scalarField mfp(nSteps_, 0.0);
    scalarField mfp2(nSteps_, 0.0);    
    scalarField time(nSteps_, 0.0);
    
    List<label> indices(nMols_, 0);
    
    scalar MFP = 0.0;
    scalar Nmols_freepaths = 0.0;
    scalar Nmols = 0.0;

    distribution velDistr(binWidthVel_);    
    distribution probDistr(binWidthProb_);

    label count = 0;
    for(label t = 0; t < nSteps_; t++)
    {
        time[t] = ( scalar(t)+1.0 )*deltaT_;
        
        if(count > 100000)
        {
            Info << "time = " << time[t] << endl;
            count = 0;
        }
        
        count++;
        
        for(label i = 0; i < nMols_; i++)
        {
            label& index = indices[i];

            scalar tStart = 0.0;             
            scalar tEnd = endTime_;
            
            if(index < endCollTimes_[i].size())
            {
                if(time[t] > endCollTimes_[i][index])
                {
                    index++;
                }
            }
            
            if(index < endCollTimes_[i].size())
            {
                tEnd = endCollTimes_[i][index];            
            }
            
            if(index < startCollTimes_[i].size()+1)
            {
                if(index != 0)
                {
                    tStart = startCollTimes_[i][index-1];
                }
            }
            
            scalar fp = 0.0; // free path
            
            if( index < collectorOfFreePaths_[i].size() )
            {
                fp = collectorOfFreePaths_[i][index];
            }
            else
            {
                fp = freePaths_[i];
            }
            
            if(time[t] > tStart)
            {
                scalar m = fp/(tEnd-tStart);
                scalar C = fp - m*tEnd;
                scalar l = m*time[t] + C;
                probDistr.add(l*rU_.refLength()/1e-9);
                MFP += l;
                Nmols_freepaths += 1.0;
                
                scalar dl = (fp*deltaT_)/(tEnd-tStart);
                scalar dV = dl/deltaT_;
                velDistr.add(dV*rU_.refVelocity());
                
//                 if(i == 0)
//                 {
//                     Info << "index = " << index 
//                         << ", tStart = " << tStart
//                         << ", tEnd = " << tEnd
//                         << ", fp = " << fp 
//                         << ", dl = " << dl
//                         << ", l = " << l
//                         << endl;
//                 }
            }
        }
        
        Nmols += scalar(nMols_);
        
        if(Nmols_freepaths > 0.0)
        {
            mfp[t] = MFP/Nmols_freepaths;
        }
        
        if(Nmols > 0.0)
        {
             mfp2[t] = MFP/Nmols;
        }
    }
    
    writeTimeData
    (
        outputPath_,
        "dropletInfo_"+fieldName_+"_method_1_cumulFreePath_nm_vs_time_ns.txt",
        time*rU_.refTime()/1e-9,
        mfp*rU_.refLength()/1e-9,
        false
    );

    writeTimeData
    (
        outputPath_,
        "dropletInfo_"+fieldName_+"_method_2_cumulFreePath_nm_vs_time_ns.txt",
        time*rU_.refTime()/1e-9,
        mfp2*rU_.refLength()/1e-9,
        false
    );    
    
    {
        
        List< Pair<scalar> > histogram = probDistr.scaledByMax();
        
        OFstream file(outputPath_/"dropletInfo_"+fieldName_+"_probability_vs_freePathDistr_nm.txt");

        if(file.good())
        {
            forAll(histogram, i)
            {
                file 
                    << histogram[i].first() << "\t"
                    << histogram[i].second() << "\t"
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void dropletInfo::write()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    {
        
        List< Pair<scalar> > histogram = velDistr.normalised();
        
        OFstream file(outputPath_/"dropletInfo_"+fieldName_+"_probabilityNormalised_vs_velocity.txt");

        if(file.good())
        {
            forAll(histogram, i)
            {
                file 
                    << histogram[i].first() << "\t"
                    << histogram[i].second() << "\t"
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void dropletInfo::write()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }    
    
}

void dropletInfo::calculateField()
{

}
    
void dropletInfo::writeField()
{
    // find the last free path that made the first start
    
    label iStart = -1;
    scalar tStart = 0.0;
    
    for (label i = 0; i < nMols_; i++)
    {
        if(startCollTimes_[i][0] > tStart)
        {
            tStart = startCollTimes_[i][0];
            iStart = i;
        }
    }
    
    Info << nl << "start at time = " << tStart << ", mol = " << iStart << endl;

    meanFreePathVsTime();
    
    // distribution for no of collisions per atom       
    {
        distribution d(binWidthColl_);
        
        for (label i = 0; i < nMols_; i++)
        {
            d.add(nCollisions_[i]);
        }        
        
        List< Pair<scalar> > histogram = d.raw();
        
        OFstream file(outputPath_/"dropletInfo_"+fieldName_+"_collisionDistribution.txt");

        if(file.good())
        {
            forAll(histogram, i)
            {
                file 
                    << histogram[i].first() << "\t"
                    << histogram[i].second() << "\t"
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void dropletInfo::write()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
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
            file >> collectorOfFreePaths_;
        }
    }
    {
        IFstream file(pathName/nameFile2_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> freePaths_;
        }
    }    
    {
        IFstream file(pathName/nameFile3_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> nCollisions_;
        }
    }
    {
        IFstream file(pathName/nameFile4_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> startCollPosX_;
        }
    }
    {
        IFstream file(pathName/nameFile5_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> startCollPosY_;
        }
    }
    {
        IFstream file(pathName/nameFile6_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> startCollPosZ_;
        }
    }
    {
        IFstream file(pathName/nameFile7_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> endCollPosX_;
        }
    }
    {
        IFstream file(pathName/nameFile8_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> endCollPosY_;
        }
    }
    {
        IFstream file(pathName/nameFile9_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> endCollPosZ_;
        }
    }    
    {
        IFstream file(pathName/nameFile10_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> startCollTimes_;
        }
    }
    {
        IFstream file(pathName/nameFile11_);

        bool goodFile = file.good();

        if(goodFile)
        {
            file >> endCollTimes_;
        }
    } 
    
//     Info << "collectorOfFreePaths = " << collectorOfFreePaths_ << endl;
//     Info << "freePaths = " << freePaths_ << endl;
//     Info << "nCollisions = " << nCollisions_ << endl;        
          
    
}

} // End namespace Foam

// ************************************************************************* //
