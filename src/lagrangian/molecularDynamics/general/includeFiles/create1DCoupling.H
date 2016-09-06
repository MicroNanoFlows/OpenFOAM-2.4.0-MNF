//
// createCoupling.H
// ~~~~~~~~~~~~
	DynamicList<couplingMD*> OneDCouplings;

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
			const PtrList<entry> cplInfo
			(
				couplingDict.lookup("couplingConfigurations")
			);

			if (cplInfo.size() == 3)
			{
				const entry& coupling = cplInfo[0];

				if (!coupling.isDict())
				{
					FatalIOErrorIn("mdFOAM", couplingDict)
									<< "Entry " << coupling << " in couplingConfigurations section is not a"
									<< " valid dictionary." << exit(FatalIOError);
				}

				word couplingName = coupling.keyword();

				if(couplingName.compare("OneD") == 0)
				{
					wordList couplingNames;
					coupling.dict().lookup("names") >> couplingNames;

					for(int j=0; j<couplingNames.size(); ++j)
					{
						Info << couplingNames[j] << endl;
						//couplingMD mdCpl1D(1, cplName);
					}
				}
				else
				{
					FatalIOErrorIn("mdFOAM", couplingDict)
									<< "Entry " << coupling << " in couplingConfigurations section is not"
									<< " named correctly (1D)." << exit(FatalIOError);
				}
			}
			else
			{
				FatalIOErrorIn("mdFOAM", couplingDict)
				                << "Incorrect number of coupling configurations, found "
								<< cplInfo.size() << ", need 3" << exit(FatalIOError);
			}
		}

		couplingDict.close();
	}
