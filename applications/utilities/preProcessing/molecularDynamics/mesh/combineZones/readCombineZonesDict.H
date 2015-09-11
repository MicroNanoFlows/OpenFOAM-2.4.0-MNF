dictionary combineRegionsDict =
    IOdictionary
    (
        IOobject
        (
            "combineZonesDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


PtrList<entry> regionsToCombineList(combineRegionsDict.lookup("zonesToCombine"));

List<bool> writeCellSets(regionsToCombineList.size(), false);
List<word> addedRegionsNames(regionsToCombineList.size());
List<bool> writeSubtractedRegions(regionsToCombineList.size(), false);
List<word> subtractedRegionsNames(regionsToCombineList.size());
List< List<word> > regionZoneNames(regionsToCombineList.size());

const cellZoneMesh& cellZones = mesh.cellZones();

forAll(regionsToCombineList, rC)
{
    const entry& regionToCombine = regionsToCombineList[rC];
    const dictionary& combineRegionDict = regionToCombine.dict();

    PtrList<entry> zoneList(combineRegionDict.lookup("zones"));
    List<word> zoneNames(zoneList.size());

    forAll(zoneList, z)
    {
        const entry& zone = zoneList[z];
        const dictionary& zoneDict = zone.dict();
    
        const word zoneName(zoneDict.lookup("name"));
    
        if(findIndex(zoneNames, zoneName) != -1)
        {
            FatalErrorIn(args.executable())
                << "Different zone name must be used: " 
                << zoneName << nl << "in: "
                << args.caseName()/runTime.system()/"combineZonesDict"
                << exit(FatalError);
        }
    
        const label& regionId = cellZones.findZoneID(zoneName);
    
        if(regionId == -1)
        {
            FatalErrorIn(args.executable())
                << "Cannot find zone: " << zoneName << " in cellZones."
                << nl << " Check: "
                << args.caseName()/runTime.system()/"combineZonesDict"
                << exit(FatalError);
        }
    
        zoneNames[z] = zoneName;
    }

    regionZoneNames[rC].setSize(zoneNames.size());
    regionZoneNames[rC] = zoneNames;

    writeCellSets[rC] = Switch(combineRegionDict.lookup("writeCellSets"));
    
    const word addedRegionName(combineRegionDict.lookup("addedZoneName"));
    addedRegionsNames[rC] = addedRegionName;
    
    writeSubtractedRegions[rC] = Switch(combineRegionDict.lookup("writeSubtractedZone"));
   
    if(writeSubtractedRegions[rC])
    {
        const word subtractedRegionName(combineRegionDict.lookup("subtractedZoneName"));
        subtractedRegionsNames[rC] = subtractedRegionName;
    }
}



