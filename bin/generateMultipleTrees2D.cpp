#include <iostream>
#include "CohabitingTrees.h"
#include "ExpandTreeHelpers.h"


template<int TDim>
void templa()
{
	CircularDomainCtrl<TDim> circ_domctr(1.0, PointD<TDim>::diagonal(0.0));

	// data vectors
	double perfs = 10000;
	std::vector<unsigned int> terms {10, 10};

	CohabitingTrees<CircularDomainCtrl<TDim>, TDim> ctree(perfs, terms, circ_domctr);

	ExpandTreeHelpers::expandCohabitingTrees(ctree, 200);

	ctree.writeTreesToXML();
}

void templamaks()
{
	std::string domain_file = "Samples/shape.png";
	ImageMaskDomainCtrl<2> mask_domctr(domain_file, 128,  PointI<2>(0,0));

	// data vectors
	double perfs = 100;
	std::vector<unsigned int> terms {200, 200};

	CohabitingTrees<ImageMaskDomainCtrl<2>, 2> ctree(perfs, terms, mask_domctr);

	ExpandTreeHelpers::expandCohabitingTrees(ctree, 200);

	ctree.exportTreesDisplays();
}

int main()
{
	templa<3>();

	//templamaks();

	return 0;
}