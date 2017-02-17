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

#include "agentCount.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(agentCount, 0);
addToRunTimeSelectionTable(agentMeasurement, agentCount, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from components
agentCount::agentCount
(
    Time& t,
    const polyMesh& mesh,
    agentCloud& cloud,
    const dictionary& dict
)
:
    agentMeasurement(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    agentIds_(),
    fieldName_(propsDict_.lookup("fieldName"))
{
   
    agentIds_.clear();

    selectAgentIds ids
    (
        cloud_.cP(),
        propsDict_
    );

    agentIds_ = ids.agentIds();

    dictionary dict2 = propsDict_.subDict("region");
    
    zone_.setBox
    (
        dict2.lookup("startPoint"),
        dict2.lookup("endPoint"),
        readScalar(dict2.lookup("radius"))  
    );        
    

    
   
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

agentCount::~agentCount()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void agentCount::createField()
{
}

void agentCount::calculateField()
{
    DynamicList<label> tNs(0);
    
    scalar nAgents = 0.0;
    
    {
        IDLList<agent>::iterator mol(cloud_.begin());
        
        for (mol = cloud_.begin(); mol != cloud_.end(); ++mol)
        {
            if(findIndex(agentIds_, mol().id()) != -1)
            {
                if(zone_.contains(mol().position()))
                {
                    label id = findIndex(tNs_, mol().trackingNumber());
                    
                    if(id == -1)
                    {
                        tNs.append(mol().trackingNumber());
                        nAgents += 1.0;
                    }
                }
            }
        }
    }
    
    if(Pstream::parRun())
    {
        // add stuff here for parallel processing 
    }
    else
    {
        nAgentsCumul_ += nAgents;
        nAgentsCumulField_.append(nAgentsCumul_);
        nAgentsInstantField_.append(nAgents);
        timeField_.append(time_.timeOutputValue());
        
        forAll(tNs, i)
        {
            tNs_.append(tNs[i]);
            times_.append(time_.timeOutputValue());
        }
        
//         Info << "times = " << times_ << endl;
//         Info << "tns = " << tNs_ << endl;
    }
}

void agentCount::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            {
                scalarField timeField(timeField_.size());
                timeField.transfer(timeField_);
                
                scalarField nAgentsCumulField(nAgentsCumulField_.size());
                nAgentsCumulField.transfer(nAgentsCumulField_);
                
                scalarField nAgentsInstantField(nAgentsInstantField_.size());
                nAgentsInstantField.transfer(nAgentsInstantField_);
                
                
                writeTimeData
                (
                    casePath_,
                    "agentCount_"+fieldName_+"_cumul_N.xy",
                    timeField,
                    nAgentsCumulField,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "agentCount_"+fieldName_+"_instant_N.xy",
                    timeField,
                    nAgentsInstantField,
                    true
                );                
            
            }
            
            {
                scalarField timeField(times_.size());
                scalarField tNs(tNs_.size());

                forAll(tNs, i)
                {
                    timeField[i]=times_[i];
                    tNs[i]=scalar(tNs_[i]);
                }
                
                writeTimeData
                (
                    casePath_,
                    "agentCount_"+fieldName_+"_trackingNumbers.xy",
                    timeField,
                    tNs
                );   
            }

        }
        
        nAgentsInstantField_.clear();
        nAgentsCumulField_.clear();
        timeField_.clear();
        
        
    }    
}



void agentCount::measureDuringForceComputation
(
    agent* molI,
    agent* molJ
)
{}

void agentCount::measureDuringForceComputationSite
(
    agent* molI,
    agent* molJ,
    label sI,
    label sJ
)
{}

// const propertyField& agentCount::fields() const
// {
//     return fields_;
// }

} // End namespace Foam

// ************************************************************************* //
