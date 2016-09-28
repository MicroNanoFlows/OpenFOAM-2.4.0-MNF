//
// createCoupling.H
// ~~~~~~~~~~~~

	{
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
						FatalIOErrorIn("mdFOAM", couplingDict)
										<< "Entry " << coupling << " in couplingConfigurations section is not a"
										<< " valid dictionary." << exit(FatalIOError);
					}

					word couplingName = coupling.keyword();

					if(couplingName.compare("OneD") == 0)
					{
						if (coupling.dict().found("names"))
						{
							wordList couplingNames;
							coupling.dict().lookup("names") >> couplingNames;

							for(int j=0; j<couplingNames.size(); ++j)
							{
								couplingInterface1d newInterface;
								newInterface.couplingName = mainCouplingName;
								newInterface.interfaceName = couplingNames[j];
								newInterface.interface = new coupling1d(couplingNames[j]);
								oneDCouplings.append(newInterface);
							}
						}
					}
					else if(couplingName.compare("TwoD") == 0)
					{
						if (coupling.dict().found("names"))
						{
							wordList couplingNames;
							coupling.dict().lookup("names") >> couplingNames;

							for(int j=0; j<couplingNames.size(); ++j)
							{
								couplingInterface2d newInterface;
								newInterface.couplingName = mainCouplingName;
								newInterface.interfaceName = couplingNames[j];
								newInterface.interface = new coupling2d(couplingNames[j]);
								twoDCouplings.append(newInterface);
							}
						}
					}
					else if(couplingName.compare("ThreeD") == 0)
					{
						if (coupling.dict().found("names"))
						{
							wordList couplingNames;
							coupling.dict().lookup("names") >> couplingNames;

							for(int j=0; j<couplingNames.size(); ++j)
							{
								couplingInterface3d newInterface;
								newInterface.couplingName = mainCouplingName;
								newInterface.interfaceName = couplingNames[j];
								newInterface.interface = new coupling3d(couplingNames[j]);
								threeDCouplings.append(newInterface);
							}
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
