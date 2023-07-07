#include <iostream>
#include "CohabitingTrees.h"
#include "ExpandTreeHelpers.h"


int main()
{
	CircularDomainCtrl<2>::TPoint p_center(0, 0);
	CircularDomainCtrl<2> circ_dom(1.0, p_center);

	// data vectors
	std::vector<double> perfs {20000, 10000};
	std::vector<unsigned int> terms {100, 100};

	CohabitingTrees<CircularDomainCtrl<2>, 2> ctree(perfs, terms, circ_dom);

	ExpandTreeHelpers::expandCohabitingTrees(ctree);

	return 0;
}