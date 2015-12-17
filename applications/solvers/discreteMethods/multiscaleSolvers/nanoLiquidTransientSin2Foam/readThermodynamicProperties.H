    Info<< "Reading thermodynamicProperties\n" << endl;

    IOdictionary thermodynamicProperties
    (
        IOobject
        (
            "thermodynamicProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rho0
    (
        thermodynamicProperties.lookup("rho0")
    );

    dimensionedScalar p0
    (
        thermodynamicProperties.lookup("p0")
    );

    dimensionedScalar molMass
    (
        thermodynamicProperties.lookup("moleculeMass")
    );

    dimensionedVector extForce
    (
        thermodynamicProperties.lookup("extForce")
    );

    dimensionedScalar a0
    (
        thermodynamicProperties.lookup("a0")
    );

    dimensionedScalar a1
    (
        thermodynamicProperties.lookup("a1")
    );

    dimensionedScalar a2
    (
        thermodynamicProperties.lookup("a2")
    );

    dimensionedScalar a3
    (
        thermodynamicProperties.lookup("a3")
    );

    dimensionedScalar a4
    (
        thermodynamicProperties.lookup("a4")
    );

    dimensionedScalar a5
    (
        thermodynamicProperties.lookup("a5")
    );

    dimensionedScalar a6
    (
        thermodynamicProperties.lookup("a6")
    );

    dimensionedScalar a7
    (
        thermodynamicProperties.lookup("a7")
    );
                           
    dimensionedScalar a8
    (
        thermodynamicProperties.lookup("a8")
    );

    dimensionedScalar a9
    (
        thermodynamicProperties.lookup("a9")
    );

    dimensionedScalar a10
    (
        thermodynamicProperties.lookup("a10")
    );

    dimensionedScalar mu6
    (
        thermodynamicProperties.lookup("mu6")
    );

    dimensionedScalar mu5
    (
        thermodynamicProperties.lookup("mu5")
    ); 

    dimensionedScalar mu4
    (
        thermodynamicProperties.lookup("mu4")
    ); 


    dimensionedScalar mu3
    (
        thermodynamicProperties.lookup("mu3")
    );

    dimensionedScalar mu2
    (
        thermodynamicProperties.lookup("mu2")
    ); 

    dimensionedScalar mu1
    (
        thermodynamicProperties.lookup("mu1")
    ); 

    dimensionedScalar mu0
    (
        thermodynamicProperties.lookup("mu0")
    ); 


//    dimensionedScalar psi
//    (
//        thermodynamicProperties.lookup("psi")
//    );

    // Density offset, i.e. the constant part of the density
//    dimensionedScalar rhoO("rhoO", rho0 - psi*p0);
