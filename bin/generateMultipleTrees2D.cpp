#include <iostream>
#include "CoronaryArteryTree.h"

template <class DomCtr, int TDim>
using TTree = CoronaryArteryTree<DomCtr, TDim>;

template<class DomCtr, int TDim>
class CohabitingTrees
{
public:
	CohabitingTrees(std::size_t tree_nb, const std::vector<double> & aPerf_vec,
					const std::vector<unsigned int> & nTerm_vec, DomCtr & dom_ctr, double radius = 1.0)
	{
		if(aPerf_vec.size() != tree_nb || nTerm_vec.size() != tree_nb)
		{
			std::cout << "Data vectors to initialize the trees don't have the correct size (" << tree_nb << ") : " << std::endl
				<< "aPerf_vec.size() = " << aPerf_vec.size() << std::endl
				<< "nTerm_vec.size() = " << nTerm_vec.size() << std::endl;
		}

		for(std::size_t i = 0; i < tree_nb; i++)
		{
			myTrees.emplace_back(aPerf_vec[i], nTerm_vec[i], dom_ctr, radius);
		}
	}

private:
	std::vector< TTree<DomCtr, TDim> > myTrees;
};



int main()
{
	CircularDomainCtrl<2>::TPoint p_center(0, 0);
	CircularDomainCtrl<2> circ_dom(1.0, p_center);

	// data vectors
	std::vector<double> perfs {20000, 10000, 5000};
	std::vector<unsigned int> terms {500, 500, 500};

	CohabitingTrees<CircularDomainCtrl<2>, 2> ctree(3, perfs, terms, circ_dom);

	return 0;
}