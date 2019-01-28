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
        bool oneDCreated = false;
        bool twoDCreated = false;
        bool threeDCreated = false;

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

            oneDInterfaces.domainName = mainCouplingName;
            twoDInterfaces.domainName = mainCouplingName;
            threeDInterfaces.domainName = mainCouplingName;

            const PtrList<entry> couplingConfigurationsList
            (
                couplingDict.lookup("couplingConfigurations")
            );

            //Iterate through configuration lists (1D/2D/3D)
            forAll(couplingConfigurationsList, cC)
            {
                const entry& couplingConfigI = couplingConfigurationsList[cC];
                const dictionary& couplingConfigIDict = couplingConfigI.dict();

                wordList iFaceTOC(couplingConfigIDict.toc());

                List<word> interfaceNames(iFaceTOC.size());
                List<List<word> > zoneNames(iFaceTOC.size());
                List<vector> domainStarts(iFaceTOC.size());
                List<vector> domainEnds(iFaceTOC.size());
                List<bool> send(iFaceTOC.size());
                List<bool> receive(iFaceTOC.size());
                List<bool> smart_send(iFaceTOC.size());

                bool domainStartFound = false;
                bool domainEndFound = false;

                //Iterate through interfaces
                forAll(iFaceTOC, cI)
                {
                    const dictionary& interfaceDict = couplingConfigIDict.subDict(iFaceTOC[cI]);

                    interfaceNames[cI] = static_cast<word>(interfaceDict.dictName());

                    if (interfaceDict.found("zones"))
                    {
                        const List<word> lclZoneNames(interfaceDict.lookup("zones"));

                        forAll(lclZoneNames, lclZones)
                        {
                            zoneNames[cI].append(lclZoneNames[lclZones]);
                        }
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "Zone missing" << exit(FatalIOError);
                    }

                    domainStarts[cI][0] = -VSMALL;
                    domainStarts[cI][1] = -VSMALL;
                    domainStarts[cI][2] = -VSMALL;
                    domainEnds[cI][0] = -VSMALL;
                    domainEnds[cI][1] = -VSMALL;
                    domainEnds[cI][2] = -VSMALL;

                    if (interfaceDict.found("domainStart"))
                    {
                        domainStartFound = true;
                        scalarList domainStartList;
                        interfaceDict.lookup("domainStart") >> domainStartList;

                        if(domainStartList.size() == 3)
                        {
                            domainStarts[cI][0] = domainStartList[0];
                            domainStarts[cI][1] = domainStartList[1];
                            domainStarts[cI][2] = domainStartList[2];
                        }
                        else
                        {
                            FatalIOErrorIn("", couplingDict)
                                           << "domainStart must contain 3 entries" << exit(FatalIOError);
                        }
                    }

                    if (interfaceDict.found("domainEnd"))
                    {
                        domainEndFound = true;
                        scalarList domainEndList;
                        interfaceDict.lookup("domainEnd") >> domainEndList;

                        if(domainEndList.size() == 3)
                        {
                            domainEnds[cI][0] = domainEndList[0];
                            domainEnds[cI][1] = domainEndList[1];
                            domainEnds[cI][2] = domainEndList[2];
                        }
                        else
                        {
                            FatalIOErrorIn("", couplingDict)
                                           << "domainEnd must contain 3 entries" << exit(FatalIOError);
                        }
                    }

                    if (interfaceDict.found("sending"))
                    {
                        send[cI] = Switch(interfaceDict.lookup("sending"));
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "Missing interface sending entry" << exit(FatalIOError);
                    }

                    if (interfaceDict.found("receiving"))
                    {
                        receive[cI] = Switch(interfaceDict.lookup("receiving"));
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "Missing interface receiving entry" << exit(FatalIOError);
                    }

                    if (interfaceDict.found("smart_send"))
                    {
                        smart_send[cI] = Switch(interfaceDict.lookup("smart_send"));
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "Missing interface smart_send entry" << exit(FatalIOError);
                    }
                }

                if(iFaceTOC.size() > 0)
                {
                    if(couplingConfigIDict.dictName() == "OneDInterfaces")
                    {
                        //Explicit domain extents defined
                        if(domainStartFound && domainEndFound)
                        {
                            oneDInterfaces.interfaces = new coupling1d(mainCouplingName, zoneNames, interfaceNames, domainStarts, domainEnds, send, receive, smart_send);
                        }
                        else //Domain extents extracted from zone
                        {
                            oneDInterfaces.interfaces = new coupling1d(mainCouplingName, zoneNames, interfaceNames, send, receive, smart_send);
                        }

                        oneDCreated = true;
                    }

                    if(couplingConfigIDict.dictName() == "TwoDInterfaces")
                    {
                        //Explicit domain extents defined
                        if(domainStartFound && domainEndFound)
                        {
                            twoDInterfaces.interfaces = new coupling2d(mainCouplingName, zoneNames, interfaceNames, domainStarts, domainEnds, send, receive, smart_send);
                        }
                        else //Domain extents extracted from zone
                        {
                            twoDInterfaces.interfaces = new coupling2d(mainCouplingName, zoneNames, interfaceNames, send, receive, smart_send);
                        }

                        twoDCreated = true;
                    }

                    if(couplingConfigIDict.dictName() == "ThreeDInterfaces")
                    {
                        //Explicit domain extents defined
                        if(domainStartFound && domainEndFound)
                        {
                            threeDInterfaces.interfaces = new coupling3d(mainCouplingName, zoneNames, interfaceNames, domainStarts, domainEnds, send, receive, smart_send);
                        }
                        else //Domain extents extracted from zone
                        {
                            threeDInterfaces.interfaces = new coupling3d(mainCouplingName, zoneNames, interfaceNames, send, receive, smart_send);
                        }

                        threeDCreated = true;
                    }
                }
            }
        }

        couplingDict.close();

        if(oneDCreated)
        {
            //Announce send and receive spans for 1D interfaces that have smart send enabled
            for(size_t i=0; i<oneDInterfaces.interfaces->size(); i++)
            {
                if(oneDInterfaces.interfaces->getInterfaceSmartSendStatus(i)) //Check if smart send should be enabled for this interface
                {
                    Foam::vector domainStart(Foam::vector::zero);
                    Foam::vector domainEnd(Foam::vector::zero);

                    if(oneDInterfaces.interfaces->getInterfaceExtentsStatus(i)) //Using zone for extents details
                    {
                        label regionId = -1;
                        const cellZoneMesh& cellZones = mesh.cellZones();
                        point min(VGREAT, VGREAT, VGREAT);
                        point max(-VSMALL, -VSMALL, -VSMALL);

                        List<word> interfaceZones = oneDInterfaces.interfaces->getInterfaceZoneNames(i);

                        forAll(interfaceZones, ifZone)
                        {
                            regionId = cellZones.findZoneID(interfaceZones[ifZone]);
                            if(regionId == -1)
                            {
                                FatalErrorIn("createCouplings.H")
                                    << "Cannot find zone: " << interfaceZones[ifZone] << nl
                                    << exit(FatalError);
                            }

                            //List of points that make up the coupling zone
                            labelList controlZoneCells = mesh.cellZones()[regionId];

                            //Iterate through coupling zone points and find overall extents
                            forAll(controlZoneCells, p)
                            {
                                const labelList& cellPts = mesh.cellPoints()[controlZoneCells[p]];

                                forAll(cellPts, pI)
                                {
                                  const point& currPoint = mesh.points()[cellPts[pI]];

                                  if(currPoint[0] < min[0]) min[0] = currPoint[0];
                                  if(currPoint[1] < min[1]) min[1] = currPoint[1];
                                  if(currPoint[2] < min[2]) min[2] = currPoint[2];
                                  if(currPoint[0] > max[0]) max[0] = currPoint[0];
                                  if(currPoint[1] > max[1]) max[1] = currPoint[1];
                                  if(currPoint[2] > max[2]) max[2] = currPoint[2];
                                }
                            }
                        }

                        domainStart[0] = min[0] - min[0] * 0.01;
                        domainStart[1] = min[1] - min[1] * 0.01;
                        domainStart[2] = min[2] - min[2] * 0.01;
                        domainEnd[0] = max[0] + max[0] * 0.01;
                        domainEnd[1] = max[1] + max[1] * 0.01;
                        domainEnd[2] = max[2] + max[2] * 0.01;
                    }
                    else //using explicit domain sizes
                    {
                        domainStart = oneDInterfaces.interfaces->getInterfaceDomainStart(i);
                        domainEnd = oneDInterfaces.interfaces->getInterfaceDomainEnd(i);
                    }

                    //Create a 1d box geometry (line)
                    mui::point1d start(domainStart[0]);
                    mui::point1d end(domainEnd[0]);
                    mui::geometry::box1d region_1d(start, end);

                    //Interface is set to send
                    if(oneDInterfaces.interfaces->getInterfaceSendStatus(i))
                    {
                        oneDInterfaces.interfaces->getInterface(i)->announce_send_span(0.0, static_cast<double>(runTime.endTime().value()), region_1d);
                    }

                    //Interface is set to receive
                    if(oneDInterfaces.interfaces->getInterfaceReceiveStatus(i))
                    {
                        oneDInterfaces.interfaces->getInterface(i)->announce_recv_span(0.0, static_cast<double>(runTime.endTime().value()), region_1d);
                    }
                }
            }
        }

        if(twoDCreated)
        {
            //Announce send and receive spans for 2D interfaces that have domain extents to enable MUI smart send
            for(size_t i=0; i<twoDInterfaces.interfaces->size(); i++)
            {
                if(twoDInterfaces.interfaces->getInterfaceSmartSendStatus(i)) //Check if smart send should be enabled for this interface
                {
                    Foam::vector domainStart(Foam::vector::zero);
                    Foam::vector domainEnd(Foam::vector::zero);

                    if(twoDInterfaces.interfaces->getInterfaceExtentsStatus(i)) //Using zone for extents details
                    {
                        label regionId = -1;
                        const cellZoneMesh& cellZones = mesh.cellZones();
                        point min(VGREAT, VGREAT, VGREAT);
                        point max(-VSMALL, -VSMALL, -VSMALL);

                        List<word> interfaceZones = twoDInterfaces.interfaces->getInterfaceZoneNames(i);

                        forAll(interfaceZones, ifZone)
                        {
                            regionId = cellZones.findZoneID(interfaceZones[ifZone]);
                            if(regionId == -1)
                            {
                                FatalErrorIn("createCouplings.H")
                                    << "Cannot find zone: " << interfaceZones[ifZone] << nl
                                    << exit(FatalError);
                            }

                            //List of points that make up the coupling zone
                            labelList controlZoneCells = mesh.cellZones()[regionId];

                            //Iterate through coupling zone points and find overall extents
                            forAll(controlZoneCells, p)
                            {
                                const labelList& cellPts = mesh.cellPoints()[controlZoneCells[p]];

                                forAll(cellPts, pI)
                                {
                                  const point& currPoint = mesh.points()[cellPts[pI]];

                                  if(currPoint[0] < min[0]) min[0] = currPoint[0];
                                  if(currPoint[1] < min[1]) min[1] = currPoint[1];
                                  if(currPoint[2] < min[2]) min[2] = currPoint[2];
                                  if(currPoint[0] > max[0]) max[0] = currPoint[0];
                                  if(currPoint[1] > max[1]) max[1] = currPoint[1];
                                  if(currPoint[2] > max[2]) max[2] = currPoint[2];
                                }
                            }
                        }

                        domainStart[0] = min[0] - min[0] * 0.01;
                        domainStart[1] = min[1] - min[1] * 0.01;
                        domainStart[2] = min[2] - min[2] * 0.01;
                        domainEnd[0] = max[0] + max[0] * 0.01;
                        domainEnd[1] = max[1] + max[1] * 0.01;
                        domainEnd[2] = max[2] + max[2] * 0.01;
                    }
                    else //using explicit domain sizes
                    {
                        domainStart = twoDInterfaces.interfaces->getInterfaceDomainStart(i);
                        domainEnd = twoDInterfaces.interfaces->getInterfaceDomainEnd(i);
                    }

                    //Create a 2d box geometry
                    mui::point2d start(domainStart[0], domainStart[1]);
                    mui::point2d end(domainEnd[0], domainEnd[1]);
                    mui::geometry::box2d region_2d(start, end);

                    //Interface is set to send
                    if(twoDInterfaces.interfaces->getInterfaceSendStatus(i))
                    {
                        twoDInterfaces.interfaces->getInterface(i)->announce_send_span(0.0, static_cast<double>(runTime.endTime().value()), region_2d);
                    }

                    //Interface is set to receive
                    if(twoDInterfaces.interfaces->getInterfaceReceiveStatus(i))
                    {
                        twoDInterfaces.interfaces->getInterface(i)->announce_recv_span(0.0, static_cast<double>(runTime.endTime().value()), region_2d);
                    }
                }
            }
        }

        if(threeDCreated)
        {
            //Announce send and receive spans for 3D interfaces that have domain extents to enable MUI smart send
            for(size_t i=0; i<threeDInterfaces.interfaces->size(); i++)
            {
                if(threeDInterfaces.interfaces->getInterfaceSmartSendStatus(i)) //Check if smart send should be enabled for this interface
                {
                    Foam::vector domainStart(Foam::vector::zero);
                    Foam::vector domainEnd(Foam::vector::zero);

                    if(threeDInterfaces.interfaces->getInterfaceExtentsStatus(i)) //Using zone for extents details
                    {
                        label regionId = -1;
                        const cellZoneMesh& cellZones = mesh.cellZones();
                        point min(VGREAT, VGREAT, VGREAT);
                        point max(-VSMALL, -VSMALL, -VSMALL);

                        List<word> interfaceZones = threeDInterfaces.interfaces->getInterfaceZoneNames(i);

                        forAll(interfaceZones, ifZone)
                        {
                            regionId = cellZones.findZoneID(interfaceZones[ifZone]);
                            if(regionId == -1)
                            {
                                FatalErrorIn("createCouplings.H")
                                    << "Cannot find zone: " << interfaceZones[ifZone] << nl
                                    << exit(FatalError);
                            }

                            //List of points that make up the coupling zone
                            labelList controlZoneCells = mesh.cellZones()[regionId];

                            //Iterate through coupling zone points and find overall extents
                            forAll(controlZoneCells, p)
                            {
                                const labelList& cellPts = mesh.cellPoints()[controlZoneCells[p]];

                                forAll(cellPts, pI)
                                {
                                  const point& currPoint = mesh.points()[cellPts[pI]];

                                  if(currPoint[0] < min[0]) min[0] = currPoint[0];
                                  if(currPoint[1] < min[1]) min[1] = currPoint[1];
                                  if(currPoint[2] < min[2]) min[2] = currPoint[2];
                                  if(currPoint[0] > max[0]) max[0] = currPoint[0];
                                  if(currPoint[1] > max[1]) max[1] = currPoint[1];
                                  if(currPoint[2] > max[2]) max[2] = currPoint[2];
                                }
                            }
                        }

                        domainStart[0] = min[0] - min[0] * 0.01;
                        domainStart[1] = min[1] - min[1] * 0.01;
                        domainStart[2] = min[2] - min[2] * 0.01;
                        domainEnd[0] = max[0] + max[0] * 0.01;
                        domainEnd[1] = max[1] + max[1] * 0.01;
                        domainEnd[2] = max[2] + max[2] * 0.01;
                    }
                    else //using explicit domain sizes
                    {
                        domainStart = threeDInterfaces.interfaces->getInterfaceDomainStart(i);
                        domainEnd = threeDInterfaces.interfaces->getInterfaceDomainEnd(i);
                    }

                    //Create a 3d box geometry
                    mui::point3d start(domainStart[0], domainStart[1], domainStart[2]);
                    mui::point3d end(domainEnd[0], domainEnd[1], domainEnd[2]);
                    mui::geometry::box3d region_3d(start, end);

                    //Interface is set to send
                    if(threeDInterfaces.interfaces->getInterfaceSendStatus(i))
                    {
                        threeDInterfaces.interfaces->getInterface(i)->announce_send_span(0.0, static_cast<double>(runTime.endTime().value()), region_3d);
                    }

                    //Interface is set to receive
                    if(threeDInterfaces.interfaces->getInterfaceReceiveStatus(i))
                    {
                        threeDInterfaces.interfaces->getInterface(i)->announce_recv_span(0.0, static_cast<double>(runTime.endTime().value()), region_3d);
                    }
                }
            }
        }
    }
