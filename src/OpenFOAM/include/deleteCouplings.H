//
// deleteCouplings.H
// ~~~~~~~~~~~~
	for(size_t i=0; i<oneDCouplings.size(); ++i)
	{
		delete oneDCouplings[i].interface;
	}

	for(size_t i=0; i<twoDCouplings.size(); ++i)
	{
		delete twoDCouplings[i].interface;
	}

	for(size_t i=0; i<threeDCouplings.size(); ++i)
	{
		delete threeDCouplings[i].interface;
	}
