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

            //Extract 1D, 2D and 3D sub dictionaries
            dictionary oneDsubDict
            (
                couplingDict.subDict("OneDInterfaces")
            );

            dictionary twoDsubDict
            (
                couplingDict.subDict("TwoDInterfaces")
            );

            dictionary threeDsubDict
            (
                couplingDict.subDict("ThreeDInterfaces")
            );

            const PtrList<dictionary> cplInfo1d
            (
                oneDsubDict.lookup("interface")
            );

            const PtrList<dictionary> cplInfo2d
            (
                twoDsubDict.lookup("interface")
            );

            const PtrList<dictionary> cplInfo3d
            (
                threeDsubDict.lookup("interface")
            );

            List<coupling1d::interface> lclOneDinterfaces(cplInfo1d.size());
            //Read 1D interfaces
            for (int i=0; i<cplInfo1d.size(); ++i)
            {
                coupling1d::interface newInterface;

                if (cplInfo1d[i].found("name"))
                {
                    newInterface.name = static_cast<word>(cplInfo1d[i].lookup("name"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Interface name missing" << exit(FatalIOError);
                }

                newInterface.domainStart[0] = -VSMALL;
                newInterface.domainStart[1] = -VSMALL;
                newInterface.domainStart[2] = -VSMALL;
                newInterface.domainEnd[0] = -VSMALL;
                newInterface.domainEnd[1] = -VSMALL;
                newInterface.domainEnd[2] = -VSMALL;

                if (cplInfo1d[i].found("domainStart"))
                {
                    scalarList domainStartList;
                    cplInfo1d[i].lookup("domainStart") >> domainStartList;

                    if(domainStartList.size() == 3)
                    {
                        newInterface.domainStart[0] = domainStartList[0];
                        newInterface.domainStart[1] = domainStartList[1];
                        newInterface.domainStart[2] = domainStartList[2];
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "domainStart must contain 3 entries" << exit(FatalIOError);
                    }
                }

                if (cplInfo1d[i].found("domainEnd"))
                {
                    scalarList domainEndList;
                    cplInfo1d[i].lookup("domainEnd") >> domainEndList;

                    if(domainEndList.size() == 3)
                    {
                        newInterface.domainEnd[0] = domainEndList[0];
                        newInterface.domainEnd[1] = domainEndList[1];
                        newInterface.domainEnd[2] = domainEndList[2];
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "domainEnd must contain 3 entries" << exit(FatalIOError);
                    }
                }

                if (cplInfo1d[i].found("sending"))
                {
                    newInterface.send = Switch(cplInfo1d[i].lookup("sending"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface sending entry" << exit(FatalIOError);
                }

                if (cplInfo1d[i].found("receiving"))
                {
                    newInterface.recieve = Switch(cplInfo1d[i].lookup("receiving"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface receiving entry" << exit(FatalIOError);
                }

                lclOneDinterfaces.append(newInterface);
            }

            oneDInterfaces.domainName = mainCouplingName;
            oneDInterfaces.interfaces = new coupling1d(mainCouplingName, lclOneDinterfaces);

            //Read 2D interfaces
            for (int i=0; i< cplInfo2d.size(); ++i)
            {
                word interfaceName;

                if (cplInfo2d[i].found("name"))
                {
                    interfaceName = static_cast<word>(cplInfo2d[i].lookup("name"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Interface name missing" << exit(FatalIOError);
                }

                vector domainStart(-VSMALL, -VSMALL, -VSMALL);
                vector domainEnd(-VSMALL, -VSMALL, -VSMALL);

                if (cplInfo2d[i].found("domainStart"))
                {
                    scalarList domainStartList;
                    cplInfo2d[i].lookup("domainStart") >> domainStartList;

                    if(domainStartList.size() == 3)
                    {
                        domainStart[0] = domainStartList[0];
                        domainStart[1] = domainStartList[1];
                        domainStart[2] = domainStartList[2];
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "domainStart must contain 3 entries" << exit(FatalIOError);
                    }
                }

                if (cplInfo2d[i].found("domainEnd"))
                {
                    scalarList domainEndList;
                    cplInfo2d[i].lookup("domainEnd") >> domainEndList;

                    if(domainEndList.size() == 3)
                    {
                        domainEnd[0] = domainEndList[0];
                        domainEnd[1] = domainEndList[1];
                        domainEnd[2] = domainEndList[2];
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "domainEnd must contain 3 entries" << exit(FatalIOError);
                    }
                }

                bool sending;

                if (cplInfo2d[i].found("sending"))
                {
                    sending = Switch(cplInfo2d[i].lookup("sending"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface sending entry" << exit(FatalIOError);
                }

                bool receiving;

                if (cplInfo2d[i].found("receiving"))
                {
                    receiving = Switch(cplInfo2d[i].lookup("receiving"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface receiving entry" << exit(FatalIOError);
                }

                //oneDInterfaces.domainName = mainCouplingName;
                //oneDInterfaces.interfaces = new coupling2d(mainCouplingName, interfaceNames);
            }

            //Read 3D interfaces
            for (int i=0; i< cplInfo3d.size(); ++i)
            {
                word interfaceName;

                if (cplInfo3d[i].found("name"))
                {
                    interfaceName = static_cast<word>(cplInfo3d[i].lookup("name"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Interface name missing" << exit(FatalIOError);
                }

                vector domainStart(-VSMALL, -VSMALL, -VSMALL);
                vector domainEnd(-VSMALL, -VSMALL, -VSMALL);

                if (cplInfo3d[i].found("domainStart"))
                {
                    scalarList domainStartList;
                    cplInfo3d[i].lookup("domainStart") >> domainStartList;

                    if(domainStartList.size() == 3)
                    {
                        domainStart[0] = domainStartList[0];
                        domainStart[1] = domainStartList[1];
                        domainStart[2] = domainStartList[2];
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "domainStart must contain 3 entries" << exit(FatalIOError);
                    }
                }

                if (cplInfo3d[i].found("domainEnd"))
                {
                    scalarList domainEndList;
                    cplInfo3d[i].lookup("domainEnd") >> domainEndList;

                    if(domainEndList.size() == 3)
                    {
                        domainEnd[0] = domainEndList[0];
                        domainEnd[1] = domainEndList[1];
                        domainEnd[2] = domainEndList[2];
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "domainEnd must contain 3 entries" << exit(FatalIOError);
                    }
                }

                bool sending;

                if (cplInfo3d[i].found("sending"))
                {
                    sending = Switch(cplInfo3d[i].lookup("sending"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface sending entry" << exit(FatalIOError);
                }

                bool receiving;

                if (cplInfo3d[i].found("receiving"))
                {
                    receiving = Switch(cplInfo3d[i].lookup("receiving"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface receiving entry" << exit(FatalIOError);
                }

                //oneDInterfaces.domainName = mainCouplingName;
                //oneDInterfaces.interfaces = new coupling2d(mainCouplingName, interfaceNames);
            }
        }

        couplingDict.close();
    }
