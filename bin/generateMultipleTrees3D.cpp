#include <iostream>

#include "CLI11.hpp"

#include "CohabitingTrees.h"
#include "ExpandTreeHelpers.h"


template<int TDim>
void templa()
{
	CircularDomainCtrl<TDim> circ_domctr(1.0, PointD<TDim>::diagonal(0.0));

	// data vectors
	double perfs = 10000;
	std::vector<unsigned int> terms {50, 50};

	CohabitingTrees<CircularDomainCtrl<TDim>, TDim> ctree(perfs, terms, circ_domctr);

	ExpandTreeHelpers::expandCohabitingTrees(ctree, 200);

	ctree.writeTreesToXML();
}

void templamaks()
{
	std::string domain_file = "Samples/shape.png";
	CircularDomainCtrl<2> circ_domctr(1.0, PointD<2>::diagonal(0.0));

	// data vectors
	double perfs = 100;
	std::vector<unsigned int> terms {4, 4};

	CohabitingTrees<CircularDomainCtrl<2>, 2> ctree(perfs, terms, circ_domctr);

	ExpandTreeHelpers::expandCohabitingTrees(ctree, 200);

	ctree.exportTreesDisplays();
}

int main()
{
	templa<3>();

	//templamaks();

	return 0;
}