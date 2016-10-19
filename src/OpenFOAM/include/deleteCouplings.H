//
// deleteCouplings.H
// ~~~~~~~~~~~~
	if(oneDInterfaces.interfaces != NULL)
	{
		delete oneDInterfaces.interfaces;
	}

	if(twoDInterfaces.interfaces != NULL)
	{
		delete twoDInterfaces.interfaces;
	}

	if(threeDInterfaces.interfaces != NULL)
	{
		delete threeDInterfaces.interfaces;
	}

	//- If this is not a parallel run then need to finalise MPI (otherwise this is handled by MUI due to the use of split_by_app() function)
	if (!args.parRunControl().parRun())
	{
		MPI_Finalize();
	}
