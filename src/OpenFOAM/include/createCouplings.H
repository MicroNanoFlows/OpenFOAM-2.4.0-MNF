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

            List<word> interfaceNames(cplInfo1d.size());
            List<word> zoneNames(cplInfo1d.size());
            List<vector> domainStarts(cplInfo1d.size());
            List<vector> domainEnds(cplInfo1d.size());
            List<bool> send(cplInfo1d.size());
            List<bool> receive(cplInfo1d.size());
            List<bool> smart_send(cplInfo1d.size());

            bool domainStartFound = false;
            bool domainEndFound = false;

            //Read 1D interfaces
            for (int i=0; i<cplInfo1d.size(); i++)
            {
                if (cplInfo1d[i].found("name"))
                {
                    interfaceNames[i] = static_cast<word>(cplInfo1d[i].lookup("name"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Interface name missing" << exit(FatalIOError);
                }

                if (cplInfo1d[i].found("zone"))
                {
                    zoneNames[i] = static_cast<word>(cplInfo1d[i].lookup("zone"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Zone missing" << exit(FatalIOError);
                }

                domainStarts[i][0] = -VSMALL;
                domainStarts[i][1] = -VSMALL;
                domainStarts[i][2] = -VSMALL;
                domainEnds[i][0] = -VSMALL;
                domainEnds[i][1] = -VSMALL;
                domainEnds[i][2] = -VSMALL;

                if (cplInfo1d[i].found("domainStart"))
                {
                    domainStartFound = true;
                    scalarList domainStartList;
                    cplInfo1d[i].lookup("domainStart") >> domainStartList;

                    if(domainStartList.size() == 3)
                    {
                        domainStarts[i][0] = domainStartList[0];
                        domainStarts[i][1] = domainStartList[1];
                        domainStarts[i][2] = domainStartList[2];
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "domainStart must contain 3 entries" << exit(FatalIOError);
                    }
                }

                if (cplInfo1d[i].found("domainEnd"))
                {
                    domainEndFound = true;
                    scalarList domainEndList;
                    cplInfo1d[i].lookup("domainEnd") >> domainEndList;

                    if(domainEndList.size() == 3)
                    {
                        domainEnds[i][0] = domainEndList[0];
                        domainEnds[i][1] = domainEndList[1];
                        domainEnds[i][2] = domainEndList[2];
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "domainEnd must contain 3 entries" << exit(FatalIOError);
                    }
                }

                if (cplInfo1d[i].found("sending"))
                {
                    send[i] = Switch(cplInfo1d[i].lookup("sending"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface sending entry" << exit(FatalIOError);
                }

                if (cplInfo1d[i].found("receiving"))
                {
                    receive[i] = Switch(cplInfo1d[i].lookup("receiving"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface receiving entry" << exit(FatalIOError);
                }

                if (cplInfo1d[i].found("smart_send"))
                {
                    smart_send[i] = Switch(cplInfo1d[i].lookup("smart_send"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface smart_send entry" << exit(FatalIOError);
                }
            }

            if(cplInfo1d.size() != 0)
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
            }

            interfaceNames.clear();
            interfaceNames.setSize(cplInfo2d.size());
            zoneNames.clear();
            zoneNames.setSize(cplInfo2d.size());
            domainStarts.clear();
            domainStarts.setSize(cplInfo2d.size());
            domainEnds.clear();
            domainEnds.setSize(cplInfo2d.size());
            send.clear();
            send.setSize(cplInfo2d.size());
            receive.clear();
            receive.setSize(cplInfo2d.size());
            smart_send.clear();
            smart_send.setSize(cplInfo2d.size());

            domainStartFound = false;
            domainEndFound = false;

            //Read 2D interfaces
            for (int i=0; i< cplInfo2d.size(); i++)
            {
                if (cplInfo2d[i].found("name"))
                {
                    interfaceNames[i] = static_cast<word>(cplInfo2d[i].lookup("name"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Interface name missing" << exit(FatalIOError);
                }

                if (cplInfo1d[i].found("zone"))
                {
                    zoneNames[i] = static_cast<word>(cplInfo1d[i].lookup("zone"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Zone missing" << exit(FatalIOError);
                }

                domainStarts[i][0] = -VSMALL;
                domainStarts[i][1] = -VSMALL;
                domainStarts[i][2] = -VSMALL;
                domainEnds[i][0] = -VSMALL;
                domainEnds[i][1] = -VSMALL;
                domainEnds[i][2] = -VSMALL;

                if (cplInfo2d[i].found("domainStart"))
                {
                    domainStartFound = true;

                    scalarList domainStartList;
                    cplInfo2d[i].lookup("domainStart") >> domainStartList;

                    if(domainStartList.size() == 3)
                    {
                        domainStarts[i][0] = domainStartList[0];
                        domainStarts[i][1] = domainStartList[1];
                        domainStarts[i][2] = domainStartList[2];
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "domainStart must contain 3 entries" << exit(FatalIOError);
                    }
                }

                if (cplInfo2d[i].found("domainEnd"))
                {
                    domainEndFound = true;

                    scalarList domainEndList;
                    cplInfo2d[i].lookup("domainEnd") >> domainEndList;

                    if(domainEndList.size() == 3)
                    {
                        domainEnds[i][0] = domainEndList[0];
                        domainEnds[i][1] = domainEndList[1];
                        domainEnds[i][2] = domainEndList[2];
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "domainEnd must contain 3 entries" << exit(FatalIOError);
                    }
                }

                if (cplInfo2d[i].found("sending"))
                {
                    send[i] = Switch(cplInfo2d[i].lookup("sending"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface sending entry" << exit(FatalIOError);
                }

                if (cplInfo2d[i].found("receiving"))
                {
                    receive[i] = Switch(cplInfo2d[i].lookup("receiving"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface receiving entry" << exit(FatalIOError);
                }

                if (cplInfo1d[i].found("smart_send"))
                {
                    smart_send[i] = Switch(cplInfo1d[i].lookup("smart_send"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface smart_send entry" << exit(FatalIOError);
                }
            }

            if(cplInfo2d.size() != 0)
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
            }

            interfaceNames.clear();
            interfaceNames.setSize(cplInfo3d.size());
            zoneNames.clear();
            zoneNames.setSize(cplInfo2d.size());
            domainStarts.clear();
            domainStarts.setSize(cplInfo3d.size());
            domainEnds.clear();
            domainEnds.setSize(cplInfo3d.size());
            send.clear();
            send.setSize(cplInfo3d.size());
            receive.clear();
            receive.setSize(cplInfo3d.size());
            smart_send.clear();
            smart_send.setSize(cplInfo2d.size());

            domainStartFound = false;
            domainEndFound = false;

            //Read 3D interfaces
            for (int i=0; i< cplInfo3d.size(); ++i)
            {
                if (cplInfo3d[i].found("name"))
                {
                    interfaceNames[i] = static_cast<word>(cplInfo3d[i].lookup("name"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Interface name missing" << exit(FatalIOError);
                }

                if (cplInfo1d[i].found("zone"))
                {
                    zoneNames[i] = static_cast<word>(cplInfo1d[i].lookup("zone"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Zone missing" << exit(FatalIOError);
                }

                domainStarts[i][0] = -VSMALL;
                domainStarts[i][1] = -VSMALL;
                domainStarts[i][2] = -VSMALL;
                domainEnds[i][0] = -VSMALL;
                domainEnds[i][1] = -VSMALL;
                domainEnds[i][2] = -VSMALL;

                if (cplInfo3d[i].found("domainStart"))
                {
                    domainStartFound = true;

                    scalarList domainStartList;
                    cplInfo3d[i].lookup("domainStart") >> domainStartList;

                    if(domainStartList.size() == 3)
                    {
                        domainStarts[i][0] = domainStartList[0];
                        domainStarts[i][1] = domainStartList[1];
                        domainStarts[i][2] = domainStartList[2];
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "domainStart must contain 3 entries" << exit(FatalIOError);
                    }
                }

                if (cplInfo3d[i].found("domainEnd"))
                {
                    domainEndFound = true;

                    scalarList domainEndList;
                    cplInfo3d[i].lookup("domainEnd") >> domainEndList;

                    if(domainEndList.size() == 3)
                    {
                        domainEnds[i][0] = domainEndList[0];
                        domainEnds[i][1] = domainEndList[1];
                        domainEnds[i][2] = domainEndList[2];
                    }
                    else
                    {
                        FatalIOErrorIn("", couplingDict)
                                       << "domainEnd must contain 3 entries" << exit(FatalIOError);
                    }
                }

                if (cplInfo3d[i].found("sending"))
                {
                    send[i] = Switch(cplInfo3d[i].lookup("sending"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface sending entry" << exit(FatalIOError);
                }

                if (cplInfo3d[i].found("receiving"))
                {
                    receive[i] = Switch(cplInfo3d[i].lookup("receiving"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface receiving entry" << exit(FatalIOError);
                }

                if (cplInfo1d[i].found("smart_send"))
                {
                    smart_send[i] = Switch(cplInfo1d[i].lookup("smart_send"));
                }
                else
                {
                    FatalIOErrorIn("", couplingDict)
                                   << "Missing interface smart_send entry" << exit(FatalIOError);
                }
            }

            if(cplInfo3d.size() != 0)
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
            }
        }

        couplingDict.close();

        //Announce send and receive spans for 1D interfaces that have smart send enabled
        for(int i=0; i<oneDInterfaces.interfaces->size(); i++)
        {
            if(oneDInterfaces.interfaces->getInterfaceSmartSendStatus(i)) //Check if smart send should be enabled for this interface
            {
                Foam::vector domainStart(Foam::vector::zero);
                Foam::vector domainEnd(Foam::vector::zero);

                if(oneDInterfaces.interfaces->getInterfaceExtentsStatus(i)) //Using zone for extents details
                {
                    //List of points that make up the coupling zone
                    labelList& controlZonePoints = mesh.pointZones()[oneDInterfaces.interfaces->getInterfaceZoneName(i)];

                    point min(VGREAT, VGREAT, VGREAT);
                    point max(-VSMALL, -VSMALL, -VSMALL);

                    //Iterate through coupling zone points and find overall extents
                    forAll(controlZonePoints, p)
                    {
                        point currVertex = mesh.points()[p];

                        if(currVertex[0] < min[0]) min[0] = currVertex[0];
                        if(currVertex[1] < min[1]) min[1] = currVertex[1];
                        if(currVertex[2] < min[2]) min[2] = currVertex[2];
                        if(currVertex[0] > max[0]) max[0] = currVertex[0];
                        if(currVertex[1] > max[1]) max[1] = currVertex[1];
                        if(currVertex[2] > max[2]) max[2] = currVertex[2];
                    }

                    domainStart[0] = min[0];
                    domainStart[1] = min[1];
                    domainStart[2] = min[2];
                    domainEnd[0] = max[0];
                    domainEnd[1] = max[1];
                    domainEnd[2] = max[2];
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

        //Announce send and receive spans for 2D interfaces that have domain extents to enable MUI smart send
        for(int i=0; i<twoDInterfaces.interfaces->size(); i++)
        {
            if(twoDInterfaces.interfaces->getInterfaceSmartSendStatus(i)) //Check if smart send should be enabled for this interface
            {
                Foam::vector domainStart(Foam::vector::zero);
                Foam::vector domainEnd(Foam::vector::zero);

                if(twoDInterfaces.interfaces->getInterfaceExtentsStatus(i)) //Using zone for extents details
                {
                    //List of points that make up the coupling zone
                    labelList& controlZonePoints = mesh.pointZones()[twoDInterfaces.interfaces->getInterfaceZoneName(i)];

                    point min(VGREAT, VGREAT, VGREAT);
                    point max(-VSMALL, -VSMALL, -VSMALL);

                    //Iterate through coupling zone points and find overall extents
                    forAll(controlZonePoints, p)
                    {
                        point currVertex = mesh.points()[p];

                        if(currVertex[0] < min[0]) min[0] = currVertex[0];
                        if(currVertex[1] < min[1]) min[1] = currVertex[1];
                        if(currVertex[2] < min[2]) min[2] = currVertex[2];
                        if(currVertex[0] > max[0]) max[0] = currVertex[0];
                        if(currVertex[1] > max[1]) max[1] = currVertex[1];
                        if(currVertex[2] > max[2]) max[2] = currVertex[2];
                    }

                    domainStart[0] = min[0];
                    domainStart[1] = min[1];
                    domainStart[2] = min[2];
                    domainEnd[0] = max[0];
                    domainEnd[1] = max[1];
                    domainEnd[2] = max[2];
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

        //Announce send and receive spans for 3D interfaces that have domain extents to enable MUI smart send
        for(int i=0; i<threeDInterfaces.interfaces->size(); i++)
        {
            if(threeDInterfaces.interfaces->getInterfaceSmartSendStatus(i)) //Check if smart send should be enabled for this interface
            {
                Foam::vector domainStart(Foam::vector::zero);
                Foam::vector domainEnd(Foam::vector::zero);

                if(threeDInterfaces.interfaces->getInterfaceExtentsStatus(i)) //Using zone for extents details
                {
                    //List of points that make up the coupling zone
                    labelList& controlZonePoints = mesh.pointZones()[threeDInterfaces.interfaces->getInterfaceZoneName(i)];

                    point min(VGREAT, VGREAT, VGREAT);
                    point max(-VSMALL, -VSMALL, -VSMALL);

                    //Iterate through coupling zone points and find overall extents
                    forAll(controlZonePoints, p)
                    {
                        point currVertex = mesh.points()[p];

                        if(currVertex[0] < min[0]) min[0] = currVertex[0];
                        if(currVertex[1] < min[1]) min[1] = currVertex[1];
                        if(currVertex[2] < min[2]) min[2] = currVertex[2];
                        if(currVertex[0] > max[0]) max[0] = currVertex[0];
                        if(currVertex[1] > max[1]) max[1] = currVertex[1];
                        if(currVertex[2] > max[2]) max[2] = currVertex[2];
                    }

                    domainStart[0] = min[0];
                    domainStart[1] = min[1];
                    domainStart[2] = min[2];
                    domainEnd[0] = max[0];
                    domainEnd[1] = max[1];
                    domainEnd[2] = max[2];
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
