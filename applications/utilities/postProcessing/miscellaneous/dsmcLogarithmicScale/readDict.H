dictionary dict =
    IOdictionary
    (
        IOobject
        (
            "logarithmicDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word inputFileName = dict.lookup("inputFileName");
    const word outputFileName = dict.lookup("outputFileName");
    
//     fileName timePath(time.path()/distributionName);

    label ROWS = readLabel(dict.lookup("rows"));
    
    ROWS--;

    scalarField time(ROWS, 0.0); 
    scalarField property(ROWS, 0.0);
    
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
                inFile >> time[r];
                inFile >> property[r];
//                 check++;
            }
        }
        
//         Info << "check on the number of rows: " << endl;
        
        
    }

    // read in list of scalar times
    List<scalar> times = List<scalar>(dict.lookup("times"));
    
    scalar deltaT = time[1]-time[0];


