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

PtrList<entry> squareList(squareZoneDict.lookup("cellZones"));

List<word> regionZoneNames(squareList.size());
List<vector> startPoints(squareList.size(), vector::zero);
List<vector> endPoints(squareList.size(), vector::zero);
List<bool> writeCellSets(squareList.size(), false);

const cellZoneMesh& cellZones = mesh.cellZones();

forAll(squareList, sL)
{
    const entry& square = squareList[sL];
    const dictionary& squareDict = square.dict();

    const word zoneName(squareDict.lookup("name"));

    if(findIndex(regionZoneNames, zoneName) != -1)
    {
        FatalErrorIn(args.executable())
            << "Different zone names must be used: " 
            << zoneName << nl << "in: "
            << args.caseName()/runTime.system()/"zonesDict"
            << exit(FatalError);
    }

    const label& regionId = cellZones.findZoneID(zoneName);

    if(regionId != -1)
    {
        FatalErrorIn(args.executable())
            << "Zone name: " << zoneName << " exists in cellZones."
            << nl << " in: "
            << args.caseName()/runTime.system()/"zonesDict"
            << nl << "Solve this by executing:"
            << nl << "> rm constant/polyMesh/*Zones"
            << nl << "> blockMesh"
            << nl << "> createCellZones"
            << exit(FatalError);
    }

    regionZoneNames[sL] = zoneName;

    startPoints[sL] = squareDict.lookup("startPoint");
    endPoints[sL] = squareDict.lookup("endPoint");

    if(squareDict.found("writeCellSets"))
    {
        writeCellSets[sL] = Switch(squareDict.lookup("writeCellSets"));        
    }
}
