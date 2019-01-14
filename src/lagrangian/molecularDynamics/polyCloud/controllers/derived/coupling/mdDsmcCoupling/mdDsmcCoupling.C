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

#include "mdDsmcCoupling.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(mdDsmcCoupling, 0);

addToRunTimeSelectionTable(polyCouplingController, mdDsmcCoupling, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
mdDsmcCoupling::mdDsmcCoupling
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict,
    couplingInterface1d &oneDInterfaces,
    couplingInterface2d &twoDInterfaces,
    couplingInterface3d &threeDInterfaces
)
:
    polyCouplingController(t, molCloud, dict, oneDInterfaces, twoDInterfaces, threeDInterfaces),
    propsDict_(dict.subDict(typeName + "Properties")),
    propsDictSend_(dict.subDict(typeName + "Sending")),
    propsDictRecv_(dict.subDict(typeName + "Receiving")),
    molIds_(),
    output_(false),
    oneDInterfaces_(),
    twoDInterfaces_(),
    threeDInterfaces_(threeDInterfaces),
#ifdef USE_MUI
    cellCentres_(),
    sendInterfaces_(),
    recvInterfaces_(),
#endif
    sending_(false),
    receiving_(false)
{
#ifdef USE_MUI
    //- Determine sending interfaces if defined
    sendInterfaces_.clear();
    sendInterfaces_.setSize(0);

    if(propsDictSend_.found("sendingInterfaces"))
    {
        const List<word> interfaces(propsDictSend_.lookup("sendingInterfaces"));
        sendInterfaces_.setSize(interfaces.size(), NULL);
        sendInterfaceNames_.setSize(interfaces.size());

        for(size_t i=0; i<interfaces.size(); ++i)
        {
            //- Find MUI interfaces
            for(size_t j=0; j<threeDInterfaces.interfaces->size(); ++j)
            {
                //- If the MUI interface is found then create a copy of its pointer address and store in sendInterfaces_
                if(threeDInterfaces.interfaces->getInterfaceName(j).compare(interfaces[i]) == 0)
                {
                    sendInterfaces_[i] = threeDInterfaces.interfaces->getInterface(j);
                    sendInterfaceNames_[i] = interfaces[i]; //- Store the receiving interface name
                    break;
                }
            }
        }

        //- Check all interfaces were found
        for(size_t i=0; i<sendInterfaces_.size(); ++i)
        {
            if(sendInterfaces_[i] == NULL)
            {
                FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
                            << "Could not find 3D MUI coupling interface (" << interfaces[i]
                            << ") to send for domain " << threeDInterfaces.domainName << exit(FatalError);
            }
            else
            {
                Info << "mdDsmcCoupling::mdDsmcCoupling(): Found 3D MUI coupling interface ("
                     << interfaces[i] << ") to send for domain " << threeDInterfaces.domainName << endl;
            }
        }
    }

    //- Determine receiving interfaces if defined
    recvInterfaces_.clear();
    recvInterfaces_.setSize(0);

    if(propsDictRecv_.found("receivingInterfaces"))
    {
        const List<word> interfaces(propsDictRecv_.lookup("receivingInterfaces"));
        recvInterfaces_.setSize(interfaces.size(), NULL);
        recvInterfaceNames_.setSize(interfaces.size());

        for(size_t i=0; i<interfaces.size(); ++i)
        {
            recvInterfaces_[i] = NULL;
            //- Find MUI interfaces
            for(size_t j=0; j<threeDInterfaces.interfaces->size(); ++j)
            {
                //- If the MUI interface is found then create a copy of its pointer address and store in sendInterfaces_
                if(threeDInterfaces.interfaces->getInterfaceName(j).compare(interfaces[i]) == 0)
                {
                    recvInterfaces_[i] = threeDInterfaces.interfaces->getInterface(j);
                    recvInterfaceNames_[i] = interfaces[i]; //- Store the receiving interface name
                    break;
                }
            }
        }

        //- Check all interfaces were found
        for(size_t i=0; i<recvInterfaces_.size(); ++i)
        {
            if(recvInterfaces_[i] == NULL)
            {
                FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
                            << "Could not find 3D MUI coupling interface (" << interfaces[i]
                            << ") to receive for domain " << threeDInterfaces.domainName << exit(FatalError);
            }
            else
            {
                Info << "mdDsmcCoupling::mdDsmcCoupling(): Found 3D MUI coupling interface ("
                     << interfaces[i] << ") to receive for domain " << threeDInterfaces.domainName << endl;
            }
        }
    }
#else
    FatalErrorIn("mdDsmcCoupling::mdDsmcCoupling()")
                << "MUI library not enabled at compilation" << exit(FatalError);
#endif

    //Determine sending properties

    if(sendInterfaces_.size() != 0)
    {
        sending_ = true;
    }

    if(recvInterfaces_.size() != 0)
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

    if(propsDict_.found("binModel"))
    {
        binModel_ =  autoPtr<binModel>
        (
            binModel::New(mesh_, propsDict_)
        );
    }

    if (propsDict_.found("output"))
    {
        output_ = Switch(propsDict_.lookup("output"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mdDsmcCoupling::~mdDsmcCoupling()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mdDsmcCoupling::initialConfiguration()
{
    sendInitialisation();
}

void mdDsmcCoupling::controlBeforeMove()
{
    receiveCoupledRegion();
}

void mdDsmcCoupling::calculateProperties()
{
    sendCoupledMolecules();
    receiveCoupledParcels();
}

void mdDsmcCoupling::sendInitialisation()
{

}

void mdDsmcCoupling::receiveCoupledRegion()
{
  /*
#ifdef USE_MUI
    //- Only receive data if at least one receiving interface is defined
    if(receiving_)
    {
        Info << threeDInterfaces_.domainName << ": Receiving MUI values for time " << time_.value()
             << " through " << recvInterfaces_.size() << " interfaces" << endl;

        if(true) //- Calculating at least one average value per cell
        {
            mui::sampler_exact3d<scalar> spatial_sampler;
            mui::chrono_sampler_exact3d chrono_sampler;

            for(size_t i=0; i<recvInterfaces_.size(); ++i) //- Iterate through the interfaces
            {
                for(int j=0; j<cellCentres_.size(); ++j) //- Iterate through the cell centres (we receive exactly as many as were sent in this example)
                {
                    if(true) //- If we are receiving mass values
                    {
                        //recvMassValues_[i][j] = recvInterfaces_[i]->fetch("m", cellCentres_[j], time_.value(), spatial_sampler, chrono_sampler);

                    }

                    if(true) //- If we are receiving density values
                    {
                        //recvDensityValues_[i][j] = recvInterfaces_[i]->fetch("p", cellCentres_[j], time_.value(), spatial_sampler, chrono_sampler);
                    }
                }
            }
        }

        //- Signal other interfaces they can move on if they have a block by committing the receive time to each interface, not needed in this example, provided for clarity
        for(size_t i=0; i<sendInterfaces_.size(); ++i)
        {
            sendInterfaces_[i]->commit(time_.value());
        }
    }
#endif
*/
}

void mdDsmcCoupling::sendCoupledMolecules()
{
  /*
#ifdef USE_MUI
    //- Only send data if at least one sending interface is defined
    if(sending_)
    {
        if(true) //- If calculating at least one average value per cell then commit initial t=0 values
        {
            scalar mass;
            scalar density;
            label molCount;
            polyMolecule* molI = NULL;
            label cellCount = 0;

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

                //- Calculate average mass for the cell
                if(molCount > 0)
                {
                    mass /= molCount;
                }

                //- Iterate through sending interfaces and push mass and/or density
                for(size_t i=0; i<sendInterfaces_.size(); ++i)
                {
                    //- Push average mass for the cell if enabled to each interface
                    if(true)
                    {
                        sendInterfaces_[i]->push("m", cellCentres_[cellCount], mass);
                    }

                    //Push cell density if enabled
                    if(true)
                    {
                        scalar volume = mesh_.V()[cellI];
                        if(volume > 0.0)
                        {
                            density = mass / mesh_.V()[cellI];
                        }
                        else
                        {
                            density = 0.0;
                        }

                        sendInterfaces_[i]->push("p", cellCentres_[cellCount], density);
                    }
                }

                cellCount++;
            }

            //Commit (transmit) values to the MUI interfaces
            for(size_t i=0; i<sendInterfaces_.size(); ++i)
            {
                sendInterfaces_[i]->commit(time_.value());
                sendInterfaces_[i]->barrier(time_.value());
            }

            Info << threeDInterfaces_.domainName << ": MUI values pushed for time " << time_.value()
                 << " to " << sendInterfaces_.size() << " interfaces" << endl;
        }
    }
#endif
*/
}

void mdDsmcCoupling::receiveCoupledParcels()
{
  /*
#ifdef USE_MUI
    //- Only receive data if at least one receiving interface is defined
    if(receiving_)
    {
        Info << threeDInterfaces_.domainName << ": Receiving MUI values for time " << time_.value()
             << " through " << recvInterfaces_.size() << " interfaces" << endl;

        if(true) //- Calculating at least one average value per cell
        {
            mui::sampler_exact3d<scalar> spatial_sampler;
            mui::chrono_sampler_exact3d chrono_sampler;

            for(size_t i=0; i<recvInterfaces_.size(); ++i) //- Iterate through the interfaces
            {
                for(int j=0; j<cellCentres_.size(); ++j) //- Iterate through the cell centres (we receive exactly as many as were sent in this example)
                {
                    if(true) //- If we are receiving mass values
                    {
                        //recvMassValues_[i][j] = recvInterfaces_[i]->fetch("m", cellCentres_[j], time_.value(), spatial_sampler, chrono_sampler);

                    }

                    if(true) //- If we are receiving density values
                    {
                        //recvDensityValues_[i][j] = recvInterfaces_[i]->fetch("p", cellCentres_[j], time_.value(), spatial_sampler, chrono_sampler);
                    }
                }
            }
        }

        //- Signal other interfaces they can move on if they have a block by committing the receive time to each interface, not needed in this example, provided for clarity
        for(size_t i=0; i<sendInterfaces_.size(); ++i)
        {
            sendInterfaces_[i]->commit(time_.value());
        }
    }
#endif
*/
}

void mdDsmcCoupling::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
  /*
#ifdef USE_MUI
    if(!Pstream::parRun() || (Pstream::parRun() && Pstream::master()))
    {
        for(int i=0; i<recvInterfaces_.size(); ++i)
        {
            fileName outputFile(timePath/recvInterfaceNames_[i]);
            OFstream of(outputFile);

            if(true) //- Output mass values
            {
                of << "Averaged mass values from coupled interface " << recvInterfaceNames_[i] << endl;
                of << "{" << endl;


                for(int j=0; j<recvMassValues_[i].size(); ++j)
                {
                   of << "(" << cellCentres_[j][0] << "," << cellCentres_[j][1] << "," << cellCentres_[j][2] << "): ";
                   of << recvMassValues_[i][j] << endl;
                }


                of << "};" << endl;
            }

            if(true) //- Output density values
            {
                of << "Averaged density values from coupled interface " << recvInterfaceNames_[i] << endl;
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
 */
}

void mdDsmcCoupling::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateCouplingControllerProperties(newDict);
    propsDict_ = newDict.subDict(typeName + "Properties");
}

} // End namespace Foam

// ************************************************************************* //
