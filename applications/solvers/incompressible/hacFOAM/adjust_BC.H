// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  adjust_BC.H
//
// Adjust boundary conditions to values from MD
//
//  Ian Cosden (cosden@seas.upenn.edu)
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



	// Get index of patch
	label inletPatchID = mesh.boundaryMesh().findPatchID("movingWall");
	
	// Get reference to boundary value
	const vectorField& faceCenters = mesh.Cf().boundaryField()[inletPatchID];
	fvPatchVectorField& hsiWallU = U.boundaryField()[inletPatchID];

//   int num_AC_bins=10;
//   scalarField vx;
	// loop over all hub faces
	forAll(hsiWallU, faceI)
	{
		// get coordinate for face center
		const vector& c = faceCenters[faceI];
		vector p;
		
		int offset;
		int out=0;
		// for(int i=0; i<num_AC_bins; i++) {
		// 	offset = i*6;
		// 	//check x location of face center
		// 	if (c[0]>AC_bin_limits[offset] && c[0]<AC_bin_limits[offset+1]) {		
		// 		//check  y location
		// 		if(c[1]>AC_bin_limits[offset+2] && c[1]<AC_bin_limits[offset+3]) {
		// 			//check z location
		// 			if(c[2]>AC_bin_limits[offset+4] && c[2]<AC_bin_limits[offset+5]) {
		// 				p = vector(vx[i]*eta_vec[0], vy[i]*eta_vec[1], vz[i]*eta_vec[2]);
		// 				hsiWallU[faceI] = p;
		// 			}
		// 		}
		// 	}  
			
		// }//end for num_AC_bins loop
		
	}//end wall face loop
