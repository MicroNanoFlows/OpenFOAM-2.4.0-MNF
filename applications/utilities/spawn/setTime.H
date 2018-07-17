//if directory exists do nothing
// else create directory and intitialse controlDict
fileName systemDir("system"); 

if( !isDir(systemDir) )
{
    mkDir(systemDir);
}


{
    OFstream file(systemDir/"controlDict");
    
    if(file.good())
    {
        file << "/*--------------------------------*- C++ -*----------------------------------*\\ " << endl;        
        file << "| =========                 |                                                 |" << endl;        
        file << "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |" << endl;        
        file << "|  \\    /   O peration     | Version:  1.6                                   |" << endl;        
        file << "|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |" << endl;        
        file << "|    \\/     M anipulation  |                                                 |" << endl;        
        file << "\\*---------------------------------------------------------------------------*/" << endl;        
        file << "FoamFile" << endl;        
        file << "{" << endl;        
        file << "     version     2.0;" << endl;        
        file << "     format      ascii;" << endl;        
        file << "     class       dictionary;" << endl;        
        file << "     object      controlDict;" << endl;        
        file << "}" << endl;        
        file << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << endl;        
        file << " " << endl;        
        file << "application     mdFoam;" << endl;        
        file << "" << endl;        
        file << "startFrom       startTime;" << endl;        
        file << "" << endl;        
        file << "startTime       0;" << endl;        
        file << "" << endl;        
        file << "stopAt          endTime;" << endl;        
        file << "" << endl;        
        file << "endTime         10;" << endl;        
        file << "" << endl;        
        file << "deltaT          0.1;" << endl;        
        file << "" << endl;        
        file << "writeControl    runTime;" << endl;        
        file << "" << endl;        
        file << "writeInterval   1;" << endl;        
        file << "" << endl;        
        file << "purgeWrite      3;" << endl;        
        file << "" << endl;        
        file << "writeFormat     ascii;" << endl;        
        file << "" << endl;        
        file << "writePrecision  12;" << endl;        
        file << "" << endl;        
        file << "writeCompression uncompressed;" << endl;        
        file << "" << endl;        
        file << "timeFormat      general;" << endl;        
        file << "" << endl;        
        file << "timePrecision   8;" << endl;        
        file << "" << endl;        
        file << "runTimeModifiable yes;" << endl;        
        file << "" << endl;        
        file << "adjustTimeStep  no;" << endl;        
        file << "" << endl;        
        file << "// ************************************************************************* //        " << endl;        

    }
}



