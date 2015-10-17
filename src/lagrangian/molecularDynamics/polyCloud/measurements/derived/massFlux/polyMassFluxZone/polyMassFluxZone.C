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

#include "polyMassFluxZone.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyMassFluxZone, 0);
addToRunTimeSelectionTable(polyField, polyMassFluxZone, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyMassFluxZone::polyMassFluxZone
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fields_(t, mesh, "dummy"),
    fieldName_(propsDict_.lookup("fieldName")),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    totalVolume_(0.0),
    length_(readScalar(propsDict_.lookup("length"))),
    unitVector_(propsDict_.lookup("unitVector")),
    mols_(0.0),
    mass_(0.0),
    mom_(vector::zero),
    velocity_(vector::zero),
    molIds_(),
	massFluxA_(1, 0.0),
	massFluxB_(1, 0.0),
	massFluxC_(1, 0.0),
	nFluxA_(1, 0.0),
    timeIndex_(0),
    nAvTimeSteps_(0.0),
    resetAtOutput_(false)
{
    resetAtOutput_ = Switch(propsDict_.lookup("resetAtOutput")); 
    
    const cellZoneMesh& cellZones = mesh_.cellZones();

    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("polyMassFluxZone::polyMassFluxZone()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

   // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    //-set the total volume
    const labelList& cells = cellZones[regionId_];

    forAll(cells, c)
    {
        const label& cellI = cells[c];
        totalVolume_ += mesh_.cellVolumes()[cellI];
    }

    if (Pstream::parRun())
    {
        reduce(totalVolume_, sumOp<scalar>());
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMassFluxZone::~polyMassFluxZone()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyMassFluxZone::createField()
{}

void polyMassFluxZone::calculateField()
{
    nAvTimeSteps_+= 1.0;
    
    const List< DynamicList<polyMolecule*> >& cellOccupancy
        = molCloud_.cellOccupancy();

    const labelList& cells = mesh_.cellZones()[regionId_];

    scalar mols = 0.0;
    scalar mass = 0.0;
    vector mom = vector::zero;


    forAll(cells, c)
    {
        const label& cellI = cells[c];
        const List<polyMolecule*>& molsInCell = cellOccupancy[cellI];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            if(findIndex(molIds_, molI->id()) != -1)
            {
                mols += 1.0;
                const scalar& massI = molCloud_.cP().mass(molI->id());
                mass += massI;
                mom += molI->v()*massI;
            }
        }
    }

    if(Pstream::parRun())
    {
        reduce(mols, sumOp<scalar>());
        reduce(mass, sumOp<scalar>());
        reduce(mom, sumOp<vector>());
    }

    if(mass > 0.0)
    {
        velocity_ += mom/mass;
    }

    mols_ += mols;
    mass_ += mass;
    mom_ += mom;


    if(time_.outputTime()) 
    {
        const scalar& nAvTimeSteps = nAvTimeSteps_;

        vector velocityA = vector::zero;
        vector velocityB = velocity_/nAvTimeSteps;

        scalar density = mass_/(length_*nAvTimeSteps);

        if(mass_ > 0)
        {
            velocityA = mom_/mass_;
        }
        
        massFluxA_[timeIndex_] = (mom_ & unitVector_)/(length_*nAvTimeSteps);

        massFluxB_[timeIndex_] = (density*velocityA & unitVector_);

        massFluxC_[timeIndex_] = (density*velocityB & unitVector_);

        nFluxA_[timeIndex_] = mols_*(velocityA & unitVector_)/(length_*nAvTimeSteps);

        if(resetAtOutput_)
        {
            mols_ = 0.0;
            mass_ = 0.0;
            mom_ = vector::zero;
            velocity_ = vector::zero;
            nAvTimeSteps_ = 1.0;
        }
    }
}

void polyMassFluxZone::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {

        if(Pstream::master())
        {
            scalarField timeField(1, runTime.timeOutputValue());


            writeTimeData
            (
                casePath_,
                "fluxZone_"+regionName_+"_"+fieldName_+"_M_A.xy",
                timeField,
                massFluxA_,
                true
            );

            writeTimeData
            (
                casePath_,
                "fluxZone_"+regionName_+"_"+fieldName_+"_M_B.xy",
                timeField,
                massFluxB_,
                true
            );

            writeTimeData
            (
                casePath_,
                "fluxZone_"+regionName_+"_"+fieldName_+"_M_C.xy",
                timeField,
                massFluxC_,
                true
            );

            writeTimeData
            (
                casePath_,
                "fluxZone_"+regionName_+"_"+fieldName_+"_N.xy",
                timeField,
                nFluxA_,
                true
            );

        	writeTimeData
			(
				casePath_,
				"fluxZone_"+regionName_+"_"+fieldName_+"_M_A.xy",
				massFluxA_,
				true
			);

			writeTimeData
			(
				casePath_,
				"fluxZone_"+regionName_+"_"+fieldName_+"_M_B.xy",
				massFluxB_,
				true
			);

			writeTimeData
			(
				casePath_,
				"fluxZone_"+regionName_+"_"+fieldName_+"_M_C.xy",
				massFluxC_,
				true
			);

			writeTimeData
			(
				casePath_,
				"fluxZone_"+regionName_+"_"+fieldName_+"_N.xy",
				nFluxA_,
				true
			);


            const reducedUnits& rU = molCloud_.redUnits();
    
            if(rU.outputSIUnits())
            {

            	writeTimeData
                (
                    casePath_,
                    "fluxZone_"+regionName_+"_"+fieldName_+"_M_A_SI.xy",
                    timeField*rU.refTime(),
                    massFluxA_*rU.refMassFlux(),
                    true
                );

                writeTimeData
                (
                    casePath_,
                    "fluxZone_"+regionName_+"_"+fieldName_+"_N_SI.xy",
                    timeField*rU.refTime(),
                    nFluxA_*rU.refMolFlux(),
                    true
                );

            	writeTimeData
				(
					casePath_,
					"fluxZone_"+regionName_+"_"+fieldName_+"_M_A_SI.xy",
					massFluxA_*rU.refMassFlux(),
					true
				);

				writeTimeData
				(
					casePath_,
					"fluxZone_"+regionName_+"_"+fieldName_+"_N_SI.xy",
					nFluxA_*rU.refMolFlux(),
					true
				);
            }
        }
    }
}

void polyMassFluxZone::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyMassFluxZone::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& polyMassFluxZone::fields() const
{
    return fields_;
}

} // End namespace Foam

// ************************************************************************* //
