/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM. 

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ioHelper.H"

#ifdef USE_MUI
    #include "mui.h"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

#ifdef USE_MUI
    if (!args.parRunControl().parRun())
    {
       MPI_Init(&argc, &argv);
    }
#endif

    Info << setprecision(10);
    reducedUnits redUnits(runTime, mesh);

    dictionary hmmDict =
    IOdictionary
    (
        IOobject
        (
            "hmmDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    ioHelper ioHelp(runTime, redUnits, hmmDict, mesh);

#ifdef USE_MUI
    Foam::word macro_domainName;
    macro_domainName="macro_0";
    std::vector<std::string> interfaceNames;
    for (label i=0; i < ioHelp.nMicro(); i++)
    {
        std::string tmp="ifs"+std::to_string(i+1);
        interfaceNames.push_back(tmp);
    }

    MPI_Comm  world = mui::mpi_split_by_app();
    std::vector<mui::uniface<mui::config_1d>* > cplMacroInterfaces=mui::create_uniface<mui::config_1d>(macro_domainName,interfaceNames);
#endif

    std::map<int,vector> a_2_c_velMap;
    std::map<int,scalarField> c_2_a_velMap_x;
    std::map<int,scalarField> c_2_a_velMap_y;
    vector vel;
    
    scalar iterationId=0;

    List<labelList> &A2C_Cells = ioHelp.cellListA2C();
    List<labelList> &C2A_Cells = ioHelp.cellListC2A();
    int nA2CDataPts_p1=A2C_Cells.size();
    int nC2ADataPts_p1=C2A_Cells.size();
    int nInterfaces=cplMacroInterfaces.size();

    Info<< "\nStarting time loop\n" << endl;
    /*
    for (size_t i=0; i < cplMacroInterfaces.size(); i++)
    {
        // send forces to LAMMPS HERE
        string lb_x="phase1_vel_x_";
	string lb_y="phase1_vel_y_";
	int lc=std::stoi(cplMacroInterfaces[i]->getIFS().substr(3,1));
	mui::point1d loc(lc);	
	double v[3];
	  
	v[0]=1.7;
	//v[1]=c_2_a_velMap_y[lc];
	
	cplMacroInterfaces[i]->push(lb_x, loc, v[0]);
	//cplMacroInterfaces[i]->push(lb_y, loc, v[1]);
	
	//Commit (transmit) values to the MUI interface			
	cplMacroInterfaces[i]->commit(iterationId);
	std::cout << "\nInitial CFD push vel " << v[0] << " ifs name: " << cplMacroInterfaces[i]->getIFS() << " location: " << lc << " iteration: " << iterationId << std::endl;
	
    }
    */

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "CourantNo.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        solve(UEqn == -fvc::grad(p));

        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            volScalarField rAU(1.0/UEqn.A());

            volVectorField HbyA("HbyA", U);
            HbyA = rAU*UEqn.H();
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                (fvc::interpolate(HbyA) & mesh.Sf())
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

	//scalar current_t = runTime.value();
        //int count = current_t/5e-13; //some tests !!!
	//if(current_t>75e-10)
	{

#ifdef USE_MUI

      #include "c_to_a.H"
      for (size_t i=0; i < cplMacroInterfaces.size(); i++)
      {
	  // send velocity to LAMMPS
	  int lc=std::stoi(cplMacroInterfaces[i]->getIFS().substr(3,1));
	  mui::point1d loc(lc);	

	  scalarField tempX=c_2_a_velMap_x[lc];
	  scalarField tempY=c_2_a_velMap_y[lc];
	  for(int j=0; j<tempX.size(); j++)
	  {
	    vel[0]=tempX[j];
	    vel[1]=tempY[j];
	    cplMacroInterfaces[i]->push("phase1_vel_x_", i*nC2ADataPts_p1+j, vel[0]);
	    cplMacroInterfaces[i]->push("phase1_vel_y_", i*nC2ADataPts_p1+j, vel[1]);
	    Info << "\nCFD push vel " << vel << " ifs name: " << cplMacroInterfaces[i]->getIFS() << " location: " << lc << " iteration: " << iterationId << endl;
	  }
	  //Commit (transmit) values to the MUI interface
	  cplMacroInterfaces[i]->commit(iterationId);
      }

      for (label i=0; i < ioHelp.nMicro(); i++)
      {
	for(int j=0;j<nA2CDataPts_p1;j++){
	  int lc=i*nA2CDataPts_p1+j;
	  vel[0]=cplMacroInterfaces[i]->fetch("phase1_md_avg_vel_x_", lc, iterationId, mui::sampler_exact1d<double>(), mui::chrono_sampler_exact1d());
	  vel[1]=cplMacroInterfaces[i]->fetch("phase1_md_avg_vel_y_", lc, iterationId, mui::sampler_exact1d<double>(), mui::chrono_sampler_exact1d());
	  vel[2]=0.0;
	  
	  a_2_c_velMap[lc] = vel;
	  Info << "\nCFD fetch vel " << vel << " ifs name: " << cplMacroInterfaces[i]->getIFS() << " location: " << lc << " iteration: " << iterationId << endl;
    
	}
      }
      
      #include "a_to_c.H"
      iterationId++;
#endif

      runTime.write();	  
      }

      Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	  << nl << endl;
    }

    Info<< "End\n" << endl;

#ifdef USE_MUI
    for(size_t i=0; i < cplMacroInterfaces.size(); i++)
    {
	delete cplMacroInterfaces[i];
    }
#endif

#ifdef USE_MUI
    if (!args.parRunControl().parRun())
    {
	MPI_Finalize();
    }
#endif


    return 0;
}


// ************************************************************************* //
