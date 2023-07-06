#include <iostream>
#include "CohabitingTrees.h"
#include "ExpandTreeHelpers.h"

class tclass
{
public:
	tclass()
	{
		for(int i = 0; i < 10; i++)
		{
			v.push_back(i);
		}
	}

	//
	void displayV()
	{
		for(int i = 0; i < 10; i++)
		{
			std::cout << v[i] << " ";
		}
		std::cout << std::endl;
	}

	std::vector<int>::iterator getB()
	{
		return v.begin();
	}

private:
	std::vector<int> v;
};

int main()
{
	/*
	CircularDomainCtrl<2>::TPoint p_center(0, 0);
	CircularDomainCtrl<2> circ_dom(1.0, p_center);

	// data vectors
	std::vector<double> perfs {20000, 10000, 5000};
	std::vector<unsigned int> terms {500, 500, 500};

	CohabitingTrees<CircularDomainCtrl<2>, 2> ctree(3, perfs, terms, circ_dom);
	*/

	tclass t;

	t.displayV();

	auto it = t.getB();
	*it = 2;

	t.displayV();

	return 0;
}