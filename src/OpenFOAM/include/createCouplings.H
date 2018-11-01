//
// createCoupling.H
// ~~~~~~~~~~~~

    {
        //- If this is not a parallel run then need to first call MPI_Init for MUI (otherwise this is called during PStream creation)
        if (!args.parRunControl().parRun())
        {
            MPI_Init(&argc, &argv);
        }

        oneDInterfaces.interfaces = NULL;
        twoDInterfaces.interfaces = NULL;
        threeDInterfaces.interfaces = NULL;

        string inputFile("couplingDict");

        IOdictionary couplingDict
        (
            IOobject
            (
                inputFile,
                runTime.system(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );

        //Read coupling dictionary file if it exists
        if (couplingDict.headerOk())
        {
            word mainCouplingName;

            if (!(couplingDict.readIfPresent("couplingName", mainCouplingName)))
            {
                FatalIOErrorIn("", couplingDict)
                                << "Missing couplingName entry" << exit(FatalIOError);
            }

            const PtrList<entry> cplInfo
            (
                couplingDict.lookup("couplingConfigurations")
            );

            if (cplInfo.size() == 3)
            {
                for (int i=0; i<cplInfo.size(); ++i)
                {
                    const entry& coupling = cplInfo[i];

                    if (!coupling.isDict())
                    {
                        FatalIOErrorIn("mdFoam+", couplingDict)
                                        << "Entry " << coupling << " in couplingConfigurations section is not a"
                                        << " valid dictionary." << exit(FatalIOError);
                    }

                    word couplingName = coupling.keyword();

                    if(couplingName.compare("OneDInterfaces") == 0)
                    {
                        if (coupling.dict().found("interfaceNames"))
                        {
                            wordList couplingNames;
                            coupling.dict().lookup("interfaceNames") >> couplingNames;

                            Foam::List<word> interfaceNames;

                            for(int j=0; j<couplingNames.size(); ++j)
                            {
                                interfaceNames.append(couplingNames[j]);
                            }

                            oneDInterfaces.domainName = mainCouplingName;
                            oneDInterfaces.interfaces = new coupling1d(mainCouplingName, interfaceNames);
                        }
                    }
                    else if(couplingName.compare("TwoDInterfaces") == 0)
                    {
                        if (coupling.dict().found("interfaceNames"))
                        {
                            wordList couplingNames;
                            coupling.dict().lookup("interfaceNames") >> couplingNames;

                            Foam::List<word> interfaceNames;

                            for(int j=0; j<couplingNames.size(); ++j)
                            {
                                interfaceNames.append(couplingNames[j]);
                            }

                            twoDInterfaces.domainName = mainCouplingName;
                            twoDInterfaces.interfaces = new coupling2d(mainCouplingName, interfaceNames);
                        }
                    }
                    else if(couplingName.compare("ThreeDInterfaces") == 0)
                    {
                        if (coupling.dict().found("interfaceNames"))
                        {
                            wordList couplingNames;
                            coupling.dict().lookup("interfaceNames") >> couplingNames;

                            Foam::List<word> interfaceNames;

                            for(int j=0; j<couplingNames.size(); ++j)
                            {
                                interfaceNames.append(couplingNames[j]);
                            }

                            threeDInterfaces.domainName = mainCouplingName;
                            threeDInterfaces.interfaces = new coupling3d(mainCouplingName, interfaceNames);
                        }
                    }
                }
            }
            else
            {
                FatalIOErrorIn("", couplingDict)
                                << "Incorrect number of coupling configurations, found "
                                << cplInfo.size() << ", need 3 (OneD, TwoD, ThreeD)" << exit(FatalIOError);
            }
        }

        couplingDict.close();
    }
