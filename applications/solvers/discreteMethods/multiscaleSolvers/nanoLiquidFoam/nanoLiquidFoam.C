/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    nanoLiquidFoam

Description
    Transient solver for trans-sonic/supersonic, laminar flow of a
    compressible liquid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvIOoptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readThermodynamicProperties.H"
	//#include "readTransportProperties.H"
    #include "createFields.H"
	#include "createFvOptions.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "compressibleCourantNo.H"

        #include "rhoEqn.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(rho, U)
          + fvm::div(phi, U)
          - fvm::laplacian(mu, U)
        );

        solve(UEqn == -fvc::grad(p));

        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            volScalarField rUA = 1.0/UEqn.A();
            U = rUA*UEqn.H();

            surfaceScalarField phid
            (
                "phid",
               fvc::interpolate(psi)
               *(
                    (fvc::interpolate(U) & mesh.Sf())
					+fvc::ddtCorr(rho, U, phi)
                  	//+ fvc::ddtCorr(rUA, rho, U, phi)
                )
            );

            phi = (fvc::interpolate(rho0 - psi*p0)/fvc::interpolate(psi))*phid;

            fvScalarMatrix pEqn
            (
                fvm::ddt(psi, p)
              + fvc::div(phi)
              + fvm::div(phid, p)
              - fvm::laplacian(rho*rUA, p)
            );

            pEqn.solve();

            phi += pEqn.flux();

            #include "compressibleContinuityErrs.H"

            U -= rUA*fvc::grad(p);
            U.correctBoundaryConditions();
        }

//        psi = (a0*pow(p,7) + b0*pow(p,6) + c0*pow(p,5) + d0*pow(p,4) + e0*pow(p,3) + f0*pow(p,2) + g0*p + h0) / (p) ;
        // psi = (a0*pow(p,3) + b0*pow(p,2) + c0*p + d0) / (p - p0) ;
 
        // rho = rho0 + psi*(p - p0);

//         psi = (a0*log(p/b0)+c0)/(p);
//	   rho =  psi*(p);
	 psi = (a0*pow(p,10) + a1*pow(p,9) + a2*pow(p,8) + a3*pow(p,7) + a4*pow(p,6) + a5*pow(p,5) 
		 + a6*pow(p,4) + a7*pow(p,3) + a8*pow(p,2)  + a9*p  + a10)/p;
	
	rho =  psi*(p);
	
	 mu = mu6 * pow(rho,6) + mu5 * pow(rho,5) + mu4 * pow(rho,4) + mu3 * pow(rho,3) +  mu2 * pow(rho,2) + mu1 * rho + mu0;
  // Info << "rho = " << rho << endl;
  // Info << "psi = " << psi << endl;

//label patchID = mesh.boundaryMesh().findPatchID("walls"); 

//const polyPatch& cPatch = mesh.boundaryMesh()[patchID];  

//vectorField nHat = cPatch.nf();
//vectorField gamma = (nHat & fvc::grad(U)().boundaryField()[patchID]) & (I - sqr(nHat));

       runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
