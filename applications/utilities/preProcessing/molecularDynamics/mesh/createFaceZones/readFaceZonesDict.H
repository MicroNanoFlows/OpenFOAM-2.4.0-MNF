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

List<bool> writeFaceSets(surfacesList.size(), false);
List<word> surfaceZoneNames(surfacesList.size());
List<vector> startingVectors(surfacesList.size(), vector::zero);
List<vector> endingVectors(surfacesList.size(), vector::zero);

const faceZoneMesh& faceZones = mesh.faceZones();


forAll(surfacesList, sL)
{
    const entry& surface = surfacesList[sL];
    const dictionary& surfaceDict = surface.dict();

    const word zoneName(surfaceDict.lookup("zoneName"));

    const label& regionId = faceZones.findZoneID(zoneName);

    if(regionId != -1)
    {
        FatalErrorIn(args.executable())
            << "FaceZone: " << zoneName << " exists on the mesh."
            << nl << " Check: "
            << args.caseName()/runTime.system()/"faceZonesDict"
            << nl << "Solve this by executing:"
            << nl << "> rm constant/polyMesh/*Zones"
            << nl << "> blockMesh"
            << nl << "> createFaceZones"            
            << exit(FatalError);
    }

    if(surfaceDict.found("writeFaceSet"))
    {
        writeFaceSets[sL] = Switch(surfaceDict.lookup("writeFaceSet"));
    }
    
    surfaceZoneNames[sL] = zoneName;

    vector startVec(surfaceDict.lookup("startPoint"));
    vector endVec(surfaceDict.lookup("endPoint"));

    startingVectors[sL] = startVec;
    endingVectors[sL] = endVec;
}





