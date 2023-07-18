#include <iostream>
#include "CohabitingTrees.h"
#include "ExpandTreeHelpers.h"

template<int TDim>
void templa()
{
	CircularDomainCtrl<TDim> circ_domctr(50, PointD<2>(0,0));

	// data vectors
	double perfs = 10000;
	std::vector<unsigned int> terms {20, 20};

	CohabitingTrees<CircularDomainCtrl<TDim>, TDim> ctree(perfs, terms, circ_domctr);

	ExpandTreeHelpers::expandCohabitingTrees(ctree, 200);

	ctree.exportTreesDisplays();
}

int main()
{
	templa<2>();

	return 0;
}