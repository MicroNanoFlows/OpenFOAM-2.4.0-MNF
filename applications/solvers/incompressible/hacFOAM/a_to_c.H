/*---------------------------------------------------------------------------*\
a_to_c.H

impose MD average velocity as boundary condition in CFD
\*---------------------------------------------------------------------------*/

// Get index of patch
label inletPatchID = mesh.boundaryMesh().findPatchID(ioHelp.OverlapInterface());
	
// Get reference to boundary value
//const vectorField& faceCenters = mesh.Cf().boundaryField()[inletPatchID];
fvPatchVectorField& hsiWallU = U.boundaryField()[inletPatchID];

vector dV;
// loop over all faces
forAll(hsiWallU, faceI)
{
  vel=a_2_c_velMap[nA2CDataPts_p1-1];  
  dV=vel-hsiWallU[faceI];
  Info << "adjust BC " << vel << " " << hsiWallU[faceI] << " " << dV  << endl;
  hsiWallU[faceI] = vector(vel[0],vel[1],0);
  //hsiWallU[faceI] += vector(dV[0],dV[1],0);
}

/*
for (int i=0; i<nInterfaces; i++)
{
  for (int j=0; j<nA2CDataPts_p1; j++)
  {
    for(int k=0; k< A2C_Cells[j].size(); k++)
    {
      int ind=A2C_Cells[j][k];

      int ij=i*nA2CDataPts_p1+j;
      U[ind].component(0)=a_2_c_velMap_x[ij];
      //Info << " TEST CFD fetch " << ind << " " << a_2_c_velMap_x[ij] << endl;
    }
  }
}
*/
