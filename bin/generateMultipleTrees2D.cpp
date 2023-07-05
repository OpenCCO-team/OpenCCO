#include <iostream>
#include "CoronaryArteryTree.h"

template<class DomCtr, int TDim>
class CohabitingTrees : public CoronaryArteryTree<DomCtr, TDim>
{
	std::vector< std::shared_ptr<CoronaryArteryTree<DomCtr, TDim> > > v;
};

int main()
{
	CircularDomainCtrl<2>::TPoint p_center(0, 0);
	CircularDomainCtrl<2> circ_dom(1.0, p_center);

	CoronaryArteryTree< CircularDomainCtrl<2>, 2 > tree1(20000, 1000, circ_dom, 1.0);
	CoronaryArteryTree< CircularDomainCtrl<2>, 2 > tree2(20000, 1000, circ_dom, 1.0);

	std::cout << tree1.myDomainController().myRadius << std::endl;

	circ_dom.myRadius = 2.0;

	std::cout << tree1.myDomainController().myRadius << std::endl;

	tree2.myDomainController().myRadius = 10000.0;

	std::cout << tree1.myDomainController().myRadius << std::endl;

	return 0;
}