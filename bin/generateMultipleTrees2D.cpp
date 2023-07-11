#include <iostream>
#include "CohabitingTrees.h"
//#include "ExpandTreeHelpers.h"


int main()
{
	ImageMaskDomainCtrl<3>::TPoint p_center(0, 0);
	ImageMaskDomainCtrl dom_ctr_mask("~/OpenCCO/Samples/tore.vol", 128, p_center, 100);

	// data vectors
	std::vector<double> perfs {20000, 10000};
	std::vector<unsigned int> terms {100, 100};

	// CohabitingTrees<CircularDomainCtrl<2>, 2> ctree(perfs, terms, circ_dom);

	// ExpandTreeHelpers::expandCohabitingTrees(ctree);

	std::vector< TPoint<3> > pts = evenlySpreadPoints<CircularDomainCtrl<3>, 3>(circ_dom,1000);

	return 0;
}