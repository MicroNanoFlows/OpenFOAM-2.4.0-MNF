dictionary planeSurfaceDict =
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


PtrList<entry> surfacesList(planeSurfaceDict.lookup("faceZones"));
List<dictionary> dictionaries(surfacesList.size());
List<word> surfaceNames(surfacesList.size());

forAll(surfacesList, sL)
{
    const entry& surface = surfacesList[sL];
    const dictionary& surfaceDict = surface.dict();
    dictionaries[sL] = surfaceDict;
    
    const word zoneName(surfaceDict.lookup("zoneName"));

    if(findIndex(surfaceNames, zoneName) != -1)
    {
        FatalErrorIn(args.executable())
            << "Different zone names must be used: " 
            << zoneName << nl << "in: "
            << args.caseName()/runTime.system()/"zonesDict"
            << exit(FatalError);
    }
    
    surfaceNames[sL] = zoneName;  
}




