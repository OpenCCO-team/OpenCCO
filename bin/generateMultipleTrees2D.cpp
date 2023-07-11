#include <iostream>
#include "CohabitingTrees.h"
#include "ExpandTreeHelpers.h"

template<int TDim>
void templa()
{
	ImageMaskDomainCtrl<TDim> dom_ctr_mask("Samples/shape.png", 128, 100);

	// data vectors
	std::vector<double> perfs {20000, 10000};
	std::vector<unsigned int> terms {100, 100};

	// CohabitingTrees<CircularDomainCtrl<TDim>, TDim> ctree(perfs, terms, circ_dom);

	// ExpandTreeHelpers::expandCohabitingTrees(ctree);

	std::vector< TPoint<TDim> > pts = evenlySpreadPoints2D(dom_ctr_mask,10);
}

int main()
{
	templa<2>();

	return 0;
}