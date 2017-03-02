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

#ifdef USE_MUI
    int nMacro=5;
	Foam::word macro_domainName;
	macro_domainName="macro_0";
	std::vector<std::string> interfaceNames;
	for(label i=0; i < nMacro; i++)
	{
		std::string tmp="ifs"+std::to_string(i+1);
		interfaceNames.push_back(tmp);
	}
	MPI_Comm  world = mui::mpi_split_by_app();
	std::vector<mui::uniface<mui::config_3d>* > cplMacroInterfaces=mui::create_uniface<mui::config_3d>(macro_domainName,interfaceNames);

#endif

    scalar nMicro=5;
    scalar iterationId=0;
    vector velocity;
    List<vector> CFD_2_MD_velocities;
    List<vector> MD_2_CFD_velocities;
	std::map<int,vector> mVelocityMap;
	double velocity_SI_2_LAMMPS=1e5;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
	    #include "adjust_BC.H"
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

#ifdef USE_MUI


		vector vel;
		for(size_t i=0; i < cplMacroInterfaces.size(); i++)
		{
			int ifsIndex=i+1;
			// convert force (SI units) to LAMMPS units (Real units)
			vel.x()=CFD_2_MD_velocities[i].x()*velocity_SI_2_LAMMPS;
			vel.y()=CFD_2_MD_velocities[i].y()*velocity_SI_2_LAMMPS;
			vel.z()=CFD_2_MD_velocities[i].z()*velocity_SI_2_LAMMPS;
			mVelocityMap[ifsIndex]=vel;
			
		}

		for(size_t i=0; i < cplMacroInterfaces.size(); i++)
		{
			// force=forces[i].x()*1.44e10;
			// send forces to LAMMPS HERE
			string lb="velocity_";
			int lc=std::stoi(cplMacroInterfaces[i]->getIFS().substr(3,1));
			mui::point3d loc;
			loc[0]=lc;
			loc[1]=lc;
			loc[2]=lc;
			double v[3];
			v[0]=mVelocityMap[lc][0];
			v[1]=mVelocityMap[lc][1];
			v[2]=mVelocityMap[lc][2];
			
			cplMacroInterfaces[i]->push(lb, loc, v[0]);
			cplMacroInterfaces[i]->push(lb, loc, v[1]);
			cplMacroInterfaces[i]->push(lb, loc, v[2]);
			//Commit (transmit) values to the MUI interface			
			cplMacroInterfaces[i]->commit(iterationId);
			std::cout << "CFD push " << "ifs : " << lc << "velocity : " << v[0]  << std::endl;

		}
  
#endif

#ifdef USE_MUI
		std::map<int,vector> velMap;

		for(label i=0; i < nMicro; i++)
		{
			int lc=std::stoi(cplMacroInterfaces[i]->getIFS().substr(3,1));
			mui::point3d loc;
			loc[0]=lc; loc[1]=lc; loc[2]=lc;
			vector v;
			v.x()=cplMacroInterfaces[i]->fetch("massrate_", loc, iterationId, mui::sampler_exact3d<double>(), mui::chrono_sampler_exact3d());
			v.y()=cplMacroInterfaces[i]->fetch("massrate_", loc, iterationId, mui::sampler_exact3d<double>(), mui::chrono_sampler_exact3d());
			v.z()=cplMacroInterfaces[i]->fetch("massrate_", loc, iterationId, mui::sampler_exact3d<double>(), mui::chrono_sampler_exact3d());
			velMap[lc] = v;
//			std::cout << "CFD fetch 0 " << "ifs : " << i << " " << lc << "mass flow rate : " << mDotMap[lc]  << " " << massrate <<std::endl;
		}

		label i=0;
		for(auto const & item : velMap) 
		{
			vector v=item.second;
			MD_2_CFD_velocities[i].x()=v.x();
			MD_2_CFD_velocities[i].y()=v.y();
			MD_2_CFD_velocities[i].z()=v.z();
			//std::cout << "CFD fetch " << "ifs : " << i << " " << item.first << " mass flow rate : " << mDot[i]  << std::endl;
			i++;
		}
		     
#endif


        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
		iterationId++;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
