    dictionary squareZoneDict =
        IOdictionary
        (
            IOobject
            (
                "zonesDict",
                runTime.system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

PtrList<entry> regionList(squareZoneDict.lookup("cellZones"));
List<dictionary> dictionaries(regionList.size());
List<word> regionZoneNames(regionList.size());

forAll(regionList, sL)
{
    const entry& region = regionList[sL];
    const dictionary& regionDict = region.dict();
    dictionaries[sL] = regionDict;
    
    const word zoneName(regionDict.lookup("zoneName"));

    if(findIndex(regionZoneNames, zoneName) != -1)
    {
        FatalErrorIn(args.executable())
            << "Different zone names must be used: " 
            << zoneName << nl << "in: "
            << args.caseName()/runTime.system()/"zonesDict"
            << exit(FatalError);
    }
    
    regionZoneNames[sL] = zoneName;
}

