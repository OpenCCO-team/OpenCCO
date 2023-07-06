#pragma once

#include <iostream>
#include "CoronaryArteryTree.h"


template <class DomCtr, int TDim>
using TTree = CoronaryArteryTree<DomCtr, TDim>;


/**
 * The goal of CohabitingTrees template class is to manage multiple trees in the same domain.
 * It can check the state of each tree, export them, etc
 * 
 * 
 */
template<class DomCtr, int TDim>
class CohabitingTrees
{
	typedef std::vector< TTree<DomCtr, TDim> > tree_vector;
	typedef typename std::vector< TTree<DomCtr, TDim> >::iterator tree_vector_iterator;

public:
	CohabitingTrees(const std::vector<double> & aPerf_vec, const std::vector<unsigned int> & nTerm_vec,
					DomCtr & dom_ctr, double radius = 1.0)
	{
		if(aPerf_vec.size() != nTerm_vec.size()
			|| aPerf_vec.empty())
		{
			std::cout << "Data vectors to construct the trees don't have the same size, or are empty : " << std::endl
				<< "\taPerf_vec.size() = " << aPerf_vec.size() << std::endl
				<< "\tnTerm_vec.size() = " << nTerm_vec.size() << std::endl;

			return;
		}

		for(std::size_t i = 0; i < aPerf_vec.size(); i++)
		{
			myTrees.emplace_back(aPerf_vec[i], nTerm_vec[i], dom_ctr, radius);
			myExpansionAttempts.push_back(0);
		}

		myIterator = myTrees.begin();
	}

	// Method
	tree_vector_iterator nextTree()
	{
		// check if myIterator is valid
		if(myIterator < myTrees.begin() || myIterator >= myTrees.end())
		{
			myIterator = myTrees.begin();
		}

		return myIterator++;
	}

	void incrementAttempt(const tree_vector_iterator & it)
	{
		std::size_t i = it - myTrees.begin();
		myExpansionAttempts[i]++;
	}

	bool expansionFinished()
	{
		bool res = true;

		for(std::size_t i = 0; i < myTrees.size(); i++)
		{
			res = res && (myExpansionAttempts[i] == myTrees[i].my_NTerm);
		}

		return res;
	}

private:
	tree_vector myTrees;
	tree_vector_iterator myIterator;
	std::vector<unsigned int> myExpansionAttempts;	// an attempt can be succesful or not; counts every try
};