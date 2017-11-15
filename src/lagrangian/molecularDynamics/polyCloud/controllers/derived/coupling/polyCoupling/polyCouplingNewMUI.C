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

#include "polyCoupling.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyCoupling, 0);

addToRunTimeSelectionTable(polyCouplingController, polyCoupling, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyCoupling::polyCoupling
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict,
    List<couplingInterface1d>& oneDCouplings,
    List<couplingInterface2d>& twoDCouplings,
    List<couplingInterface3d>& threeDCouplings
)
:
    polyCouplingController(t, molCloud, dict, oneDCouplings, twoDCouplings, threeDCouplings),
    propsDict_(dict.subDict(typeName + "Properties")),
    propsDictSend_(dict.subDict(typeName + "Sending")),
    propsDictRecv_(dict.subDict(typeName + "Receiving")),
    molIds_(),
    output_(false),
    oneDCouplings_(),
    twoDCouplings_(),
    threeDCouplings_(threeDCouplings),
    cplInterfaceName_(),
#ifdef USE_MUI
    currInterface(NULL),
    cellCentres_(),
#endif
    sending_(false),
    receiving_(false),
    sendMass_(false),
    sendDensity_(false),
    recvInterfaces_(),
    recvMass_(false),
    recvDensity_(false),
    recvMassValues_(),
    recvDensityValues_()
{
    propsDict_.readIfPresent("interfaceName", cplInterfaceName_);

    //Find MUI interface
#ifdef USE_MUI
    for(int i=0; i<threeDCouplings_.size(); ++i)
    {
        if(threeDCouplings_[i].interfaceName.compare(cplInterfaceName_) == 0)
        {
            currInterface = threeDCouplings_[i].interface->getInterface();
            break;
        }
    }

    if(currInterface == NULL)
    {
        FatalErrorIn("polyCoupling::polyCoupling()")
                    << "Could not find 3D MUI coupling interface " << cplInterfaceName_
                    << exit(FatalError);
    }
    else
    {
        Info << "Found 3D MUI coupling interface " << cplInterfaceName_ << endl;
    }
#else
    FatalErrorIn("polyCoupling::polyCoupling()")
                << "MUI library not enabled at compilation"
                << exit(FatalError);
#endif

    //Determine sending properties
    if(propsDictSend_.found("mass"))
    {
        sendMass_ = Switch(propsDictSend_.lookup("mass"));
    }

    if(propsDictSend_.found("density"))
    {
        sendDensity_ = Switch(propsDictSend_.lookup("density"));
    }

    if(sendMass_ || sendDensity_)
    {
        sending_ = true;
    }

    //Determine receiving properties
    recvInterfaces_.clear();
    const List<word> interfaces(propsDictRecv_.lookup("receivingInterfaces"));
    recvInterfaces_.setSize(interfaces.size());

    for(int i=0; i<interfaces.size(); ++i)
    {
        recvInterfaces_[i] = interfaces[i];
    }

    recvMassValues_.setSize(recvInterfaces_.size());
    recvDensityValues_.setSize(recvInterfaces_.size());

    if(propsDictRecv_.found("mass"))
    {
        recvMass_ = Switch(propsDictRecv_.lookup("mass"));
    }

    if(propsDictRecv_.found("density"))
    {
        recvDensity_ = Switch(propsDictRecv_.lookup("density"));
    }

    if(recvMass_ || recvDensity_)
    {
        receiving_ = true;
    }

    writeInTimeDir_ = true;
    writeInCase_ = true;

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    if (propsDict_.found("output"))
    {
        output_ = Switch(propsDict_.lookup("output"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyCoupling::~polyCoupling()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyCoupling::initialConfiguration()
{
#ifdef USE_MUI
    //Only send initial list of points if this interface is sending
    if(sending_)
    {
        cellCentres_.setSize(controlZone().size());
        label cellCount = 0;

        // Label list of the face labels
        const faceList &ff = mesh_.faces();
        // The coordinate sets for the individual points
        const pointField &pp = mesh_.points();

        //Start by sending list of points at cell centres in the control zone
        forAll(controlZone(), c)
        {
            const label& cellI = controlZone()[c];

            mesh_.cells()[cellI].labels(ff);
            mesh_.cells()[cellI].points(ff, pp);

            point cellCentre = mesh_.cells()[cellI].centre(pp, ff);

            //Create MUI point from OpenFOAM class
            mui::point3d centre(cellCentre.x(), cellCentre.y(), cellCentre.z());

            cellCentres_[cellCount] = centre;

            //Push the cell centre to the local MUI interface
            currInterface->push(centre);

            cellCount++;
        }

        Info << "Static point list stored to 3D MUI coupling interface "
             << cplInterfaceName_ << endl;

        for(int i=0; i<recvInterfaces_.size(); ++i)
        {
            recvMassValues_[i].setSize(cellCentres_.size());
            recvDensityValues_[i].setSize(cellCentres_.size());
        }

        if(sendMass_ || sendDensity_) //If calculating at least one average value per cell then commit initial t=0 values
        {
            scalar mass;
            scalar density;
            label molCount;
            polyMolecule* molI = NULL;
            label cellCount = 0;

            word massLbl = "m_"+cplInterfaceName_;
            word densityLbl = "p_"+cplInterfaceName_;

            forAll(controlZone(), c)
            {
                mass = 0.0;
                molCount = 0;
                const label& cellI = controlZone()[c];
                const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];

                forAll(molsInCell, m)
                {
                    molI = molsInCell[m];
                    mass += molCloud_.cP().mass(molI->id());
                    molCount++;
                }

                //Calculate average mass for the cell
                mass /= molCount;

                //Push average mass for the cell if enabled (static cell centre list used)
                if(sendMass_)
                {
                    currInterface->push(massLbl, mass);
                }

                //Push cell density if enabled (static cell centre list used)
                if(sendDensity_)
                {
                    density = mass / mesh_.V()[cellI];
                    currInterface->push(densityLbl, density);
                }

                cellCount++;
            }

            //Commit (transmit) values to the MUI interface
            currInterface->commit(time_.value());
            Info << "MUI values pushed for time " << time_.value() << endl;
            currInterface->barrier(time_.value());
        }
    }
#endif
}

void polyCoupling::sendCoupling()
{
#ifdef USE_MUI
    //Only send data if this interface is sending
    if(sending_)
    {
        if(sendMass_ || sendDensity_) //Calculating at least one average value per cell
        {
            scalar mass;
            scalar density;
            label molCount;
            polyMolecule* molI = NULL;
            label cellCount = 0;

            word massLbl = "m_"+cplInterfaceName_;
            word densityLbl = "p_"+cplInterfaceName_;

            forAll(controlZone(), c)
            {
                mass = 0.0;
                molCount = 0;
                const label& cellI = controlZone()[c];
                const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];

                forAll(molsInCell, m)
                {
                    molI = molsInCell[m];
                    mass += molCloud_.cP().mass(molI->id());
                    molCount++;
                }

                //Calculate average mass for the cell
                mass /= molCount;

                //Push average mass for the cell if enabled (static cell centre list used)
                if(sendMass_)
                {
                    currInterface->push(massLbl, mass);
                }

                //Push cell density if enabled (static cell centre list used)
                if(sendDensity_)
                {
                    density = mass / mesh_.V()[cellI];
                    currInterface->push(densityLbl, density);
                }

                cellCount++;
            }

            //Commit (transmit) values to the MUI interface
            currInterface->commit(time_.value());
            Info << "MUI values pushed for time " << time_.value() << endl;
        }
    }
#endif
}

void polyCoupling::receiveCoupling()
{
#ifdef USE_MUI
    //Only receive data if this interface is receiving
    if(receiving_)
    {
        Info << "Receiving MUI values for time " << time_.value() << endl;

        if(recvMass_ || recvDensity_) //Calculating at least one average value per cell
        {
            for(int i=0; i<recvInterfaces_.size(); ++i) //Iterate through the interfaces
            {
                word massLbl = "m_"+recvInterfaces_[i];
                word densityLbl = "p_"+recvInterfaces_[i];

                for(int j=0; j<cellCentres_.size(); ++j)
                {
                    if(recvMass_)
                    {
                        recvMassValues_[i][j] = currInterface->fetch(massLbl, cellCentres_[j], time_.value(),
                                                          mui::sampler_exact3d<scalar>(),
                                                          mui::chrono_sampler_exact3d());

                    }

                    if(recvMass_)
                    {
                        recvDensityValues_[i][j] = currInterface->fetch(densityLbl, cellCentres_[j], time_.value(),
                                                          mui::sampler_exact3d<scalar>(),
                                                          mui::chrono_sampler_exact3d());
                    }
                }
            }
            currInterface->commit(time_.value()); // Signal other interfaces they can move on if they have a block
        }
    }
#endif
}

void polyCoupling::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
#ifdef USE_MUI
    if(Pstream::master())
    {
        for(int i=0; i<recvInterfaces_.size(); ++i)
        {
            fileName outputFile(timePath/recvInterfaces_[i]);
            OFstream of(outputFile);

            if(recvMass_) //Output mass values
            {
                of << "Averaged mass values from coupled interface " << recvInterfaces_[i] << endl;
                of << "{" << endl;

                for(int j=0; j<recvMassValues_[i].size(); ++j)
                {
                    of << "(" << cellCentres_[j][0] << "," << cellCentres_[j][1] << "," << cellCentres_[j][2] << "): ";
                    of << recvMassValues_[i][j] << endl;
                }

                of << "};" << endl;
            }

            if(recvDensity_) //Output density values
            {
                of << "Averaged density values from coupled interface " << recvInterfaces_[i] << endl;
                of << "{" << endl;

                for(int j=0; j<recvDensityValues_[i].size(); ++j)
                {
                    of << "(" << cellCentres_[j][0] << "," << cellCentres_[j][1] << "," << cellCentres_[j][2] << "): ";
                    of << recvDensityValues_[i][j] << endl;
                }

                of << "};" << endl;
            }
        }
    }
#endif
}

void polyCoupling::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateCouplingControllerProperties(newDict);
    propsDict_ = newDict.subDict(typeName + "Properties");
}

} // End namespace Foam

// ************************************************************************* //
