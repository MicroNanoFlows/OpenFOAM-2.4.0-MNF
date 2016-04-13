dictionary timeCorrelationDict =
    IOdictionary
    (
        IOobject
        (
            "timeCorrelationDict",
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
    
    ROWS--;
    
    scalarField R(ROWS, 0.0);
    
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
//                 check++;
            }
            else
            {
                inFile >> temporary;
                inFile >> R[r];
//                 check++;
            }
        }
        
//         Info << "check on the number of rows: " << endl;
        
        
    }
    
    label tau = readLabel(timeCorrelationDict.lookup("tau"));
    
    label fftpoints = readLabel(timeCorrelationDict.lookup("fftpoints"));




