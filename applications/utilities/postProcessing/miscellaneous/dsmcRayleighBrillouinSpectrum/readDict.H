dictionary timeCorrelationDict =
    IOdictionary
    (
        IOobject
        (
            "rayleighBrillouinSpectrumDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word inputFileName = timeCorrelationDict.lookup("inputFileName");
    const word outputFileName = timeCorrelationDict.lookup("outputFileName");
    
//     fileName timePath(time.path()/distributionName);

    label ROWS = readLabel(timeCorrelationDict.lookup("rows"));
    
    scalar omegaMax = readScalar(timeCorrelationDict.lookup("omegaMax"));
    label omegaSteps = readLabel(timeCorrelationDict.lookup("omegaSteps"));
    
    scalar omegaJump = omegaMax/omegaSteps;
    
    ROWS--;
    
    scalarField time(ROWS, 0.0);
    scalarField Rreal(ROWS, 0.0);
    scalarField Rimaginary(ROWS, 0.0);
    
    {
        ifstream inFile;
        scalar temporary;
        
        inFile.open(inputFileName.c_str());
        
//         label check = 0;
        
        if (!inFile)
        {
            FatalErrorIn(args.executable())
                << "Unable to open file"
                << exit(FatalError);
        }

        for(label r=0; r<ROWS; r++) 
        {
            if(r == 0)
            {
                inFile >> temporary;
                inFile >> temporary;
                inFile >> temporary;
//                 check++;
            }
            else
            {
                inFile >> time[r];
                inFile >> Rreal[r];
                inFile >> Rimaginary[r];
//                 check++;
            }
        }
        
//         Info << "check on the number of rows: " << endl;
        
        
    }
    
//     label tau = readLabel(timeCorrelationDict.lookup("tau"));
    
//     label fftpoints = readLabel(timeCorrelationDict.lookup("fftpoints"));




