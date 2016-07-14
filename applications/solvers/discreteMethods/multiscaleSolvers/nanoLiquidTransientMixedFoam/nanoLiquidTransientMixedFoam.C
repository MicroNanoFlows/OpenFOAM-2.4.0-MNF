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
#include "OFstream.H"

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

//	dimensionedScalar period = dimensionedScalar("period", dimensionSet(0,0,1,0,0,0,0), 0.22e-9);
//dimensionedScalar tRed =  dimensionedScalar("tRed", dimensionSet(0,0,1,0,0,0,0), 2.16059e-12);

//% STARTING PERIOD
//dimensionedScalar P_S =  dimensionedScalar("P_S", dimensionSet(0,0,1,0,0,0,0), 100* 2.16059e-12);
scalar P_S = 100; // *2.16059e-12;
//scalar P_S = 100*2.16059e-12;
//% END PERIOD
//dimensionedScalar P_E =  dimensionedScalar("P_E", dimensionSet(0,0,1,0,0,0,0), 5000* 2.16059e-12);
scalar P_E = 5000; //*2.16059e-12;
//scalar P_E = 5000*2.16059e-12;
//dimensionedScalar endTime =  dimensionedScalar("endTime", dimensionSet(0,0,1,0,0,0,0), 10000* 2.16059e-12); 
scalar endTime = 10000; //*2.16059e-12;
//scalar endTime = 10000*2.16059e-12;
   // QUICKER TO BEGIN WITH AND HIGHER ORDER POLYNOMIAL
scalar n = 30; //*2.16059e-12;

//dimensionedScalar phase = 0; 
scalar phase = 0;
//scalar period = 0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "compressibleCourantNo.H"

        #include "rhoEqn.H"

scalar T = runTime.time().value() / 2.16059e-12;

   scalar temp = ( Foam::pow(( (T - endTime)/endTime), n));
Info << "temp = " << temp << endl;
//   dimensionedScalar freq =  dimensionedScalar("freq", dimensionSet(0,0,-1,0,0,0,0), 0); 
   scalar freq= ((1/P_E) - ( ( (1/P_E) - (1/P_S) ) * temp )) /2.16059e-12 ;
   
//   dimensionedScalar period = dimensionedScalar("period", dimensionSet(0,0,1,0,0,0,0), 1/freq);
//   dimensionedScalar period =  1/freq;
scalar  period = 1/freq;
//   phase +=  (freq * runTime.deltaT().value()/2.16059e-12)  ; // 0.0025*2.16059e-12;
Info << "phase = " << phase << endl;
   phase +=  ((freq) * runTime.deltaT().value())  ; // 0.0025*2.16059e-12;
Info << "freq = " << freq << endl;
Info << "phase = " << phase << endl;
Info << "period = " << period << endl;

   extForceField = extForce * rho / molMass * Foam::sin(2*3.1415*phase); //  *phase) ;


        fvVectorMatrix UEqn
        (
            fvm::ddt(rho, U)
          + fvm::div(phi, U)
          - fvm::laplacian(mu, U)
	  - extForceField 
        );

        solve(UEqn == -fvc::grad(p));

        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {

	 psi = (a0*pow(p,10) + a1*pow(p,9) + a2*pow(p,8) + a3*pow(p,7) + a4*pow(p,6) + a5*pow(p,5) 
		 + a6*pow(p,4) + a7*pow(p,3) + a8*pow(p,2)  + a9*p  + a10)/p;
	
	 rho =  psi*(p);
	
	 mu = mu6 * pow(rho,6) + mu5 * pow(rho,5) + mu4 * pow(rho,4) + mu3 * pow(rho,3) +  mu2 * pow(rho,2) + mu1 * rho + mu0;


            volScalarField rUA = 1.0/UEqn.A();
            U = rUA*UEqn.H();

            surfaceScalarField phid
            (
                "phid",
               fvc::interpolate(psi)
               *(
                    (fvc::interpolate(U) & mesh.Sf())
					+ fvc::ddtCorr(rho, U, phi)
                  //+ fvc::ddtPhiCorr(rUA, rho, U, phi)
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

//       rhoU =    ( rho*U.component(0) ) * (mesh.Sf()); //  fvc::snGrad(fvc::reconstruct(phi)); //   fvc::reconstruct(phi);


       runTime.write();
//	#include "binMassFluxNew.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
