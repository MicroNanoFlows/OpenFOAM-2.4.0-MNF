/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2008-2009 OpenCFD Ltd.
    \\/      M anipulation   |
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

Application
    icoMacroFOAM

Description
    Incompressible macro solver for the Internal-flow Multiscale Method

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "icoMacroModel.H"

#ifdef USE_MUI
    #include "fvCoupling.H"
//    #include "coupling1d.H"
    #include "mui.h"
#endif

using namespace Foam;

int main(int argc, char *argv[])
{
	
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

#ifdef USE_MUI
	if (!args.parRunControl().parRun())
	{
		MPI_Init(&argc, &argv);
	}
#endif

    reducedUnits redUnits(runTime, mesh);

    dictionary immDict =
    IOdictionary
    (
        IOobject
        (
            "immDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

	icoMacroModel macro(runTime, redUnits, immDict);

#ifdef USE_MUI
	Foam::word macro_domainName;
	macro_domainName="macro_0";
	std::vector<std::string> interfaceNames;
	for(label i=0; i < macro.nMicro(); i++)
	{
		std::string tmp="ifs"+std::to_string(i+1);
		interfaceNames.push_back(tmp);
		//std::cout << "cfd ifs: " << interfaceNames[i] << std::endl;
	}
	MPI_Comm  world = mui::mpi_split_by_app();
	std::vector<mui::uniface<mui::config_1d>* > cplMacroInterfaces=mui::create_uniface<mui::config_1d>(macro_domainName,interfaceNames);
	//std::cout << "TT " << interfaceNames.size() << std::endl;
#endif

    Info << "\nStarting time loop\n" << endl;

//    return 0;

	scalar force=0;
	scalar massrate=0.0;

	// set and get forces 
	List<vector> forces = macro.setAndGetForces();
/*
	for(size_t i=0; i < cplMacroInterfaces.size(); i++)
	{
		force=forces[i].x()*1.44e10; // convert to LAMMPS force units
		// send forces to LAMMPS HERE
		int lc=std::stoi(cplMacroInterfaces[i]->getIFS().substr(3,1));
		mui::point1d loc(lc);
		cplMacroInterfaces[i]->push("force_", loc, force);

		//Commit (transmit) values to the MUI interface			
		cplMacroInterfaces[i]->commit(0);
		std::cout << "CFD push 0 " << "ifs : " << lc << "force : " << force  << std::endl;
	}
*/
	//mui::point3d force;
	for(label n = 0; n < macro.nIter(); n++)
    {
        // set and get forces 
        List<vector> forces = macro.setAndGetForces();
       
#ifdef USE_MUI

		std::map<int,scalar> mForceMap;
		for(size_t i=0; i < cplMacroInterfaces.size(); i++)
		{
			int ifsIndex=i+1;
			// convert force (SI units) to LAMMPS units (Real units)
			mForceMap[ifsIndex]=forces[i].x()*1.44e10;
		}
		for(size_t i=0; i < cplMacroInterfaces.size(); i++)
		{
			// force=forces[i].x()*1.44e10;
			// send forces to LAMMPS HERE
			string lb="force_";
			int lc=std::stoi(cplMacroInterfaces[i]->getIFS().substr(3,1));
			mui::point1d loc(lc);
			force=mForceMap[lc]; 
			cplMacroInterfaces[i]->push(lb, loc, force);
//			cplMacroInterfaces[i]->push("loc_", 0, i+1);
			//Commit (transmit) values to the MUI interface			
			cplMacroInterfaces[i]->commit(n);
			std::cout << "CFD push " << "ifs : " << lc << "force : " << force  << std::endl;
		}
    
        //  RUN MD simulation HERE
#else
		macro.pseudoMD();
#endif

        List<scalar> mDot(macro.nMicro(),0.0);        
        // get mass flow rates from MD
        // and set them in icoMacroModel 
        // change from LAMMPS units to SI
#ifdef USE_MUI
		std::map<int,scalar> mDotMap;

		for(label i=0; i < macro.nMicro(); i++)
		{
			int lc=std::stoi(cplMacroInterfaces[i]->getIFS().substr(3,1));
			mui::point1d loc(lc);
			massrate = cplMacroInterfaces[i]->fetch("massrate_", loc, n, mui::sampler_exact1d<double>(), mui::chrono_sampler_exact1d());
			mDotMap[lc] = massrate*1.0;
			std::cout << "CFD fetch 0 " << "ifs : " << i << " " << lc << "mass flow rate : " << mDotMap[lc]  << " " << massrate <<std::endl;
			//mDot[i]=0.01;
		}

		//for(label i=0; i < macro.nMicro(); i++)
		label i=0;
		for(auto const & item : mDotMap) 
		{
			mDot[i]=item.second;
			std::cout << "CFD fetch " << "ifs : " << i << " " << item.first << " mass flow rate : " << mDot[i]  << std::endl;
			i++;
		}
		     
		//Commit (transmit) values to the MUI interface


        //  RUN MD simulation HERE
#endif
		macro.massFlowRates()=mDot;
        
        macro.solve();
        
        macro.write();
    }

//    Info << "End 0\n" << endl;    
#ifdef USE_MUI
//     #include "deleteCouplings.H"
	for(size_t i=0; i < cplMacroInterfaces.size(); i++)
		delete cplMacroInterfaces[i];

#endif

	//  Info << "End 1\n" << endl;

#ifdef USE_MUI
	if (!args.parRunControl().parRun())
	{
		MPI_Finalize();
	}
#endif
    return 0;
}












