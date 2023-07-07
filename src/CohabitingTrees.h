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

			// initialize first segment of trees so that they're not overlapping
			//initializeFirstSegments
		}

		myIterator = myTrees.begin();
	}

	// Method
	void nextTree()
	{
		myIterator++;

		// check if myIterator is valid
		if(myIterator < myTrees.begin() || myIterator >= myTrees.end())
		{
			myIterator = myTrees.begin();
		}
	}

	TTree<DomCtr, TDim> getCurrentTreeCopy()
	{
		return *myIterator;
	}

	void incrementAttempt()
	{
		std::size_t i = myIterator - myTrees.begin();
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

	typename TTree<DomCtr, TDim>::TPointD generateNewLocation(unsigned int nb_trials)
	{
		typename TTree<DomCtr, TDim>::TPointD res;
		double dist_threshold = myIterator->getDistanceThreshold();
		
		bool location_found = false;
		unsigned int i = 0;
		while(i < nb_trials && !location_found)
		{
			res = myIterator->myDomainController().randomPoint();

			if(validNewLocation(res, dist_threshold))
			{
				location_found = true;
			}

			// if we reach end of the loop before finding the new location,
			// loop again with less restrictive distance threshold
			i++;			
			if(!location_found && i == nb_trials)
			{
				i = 0;
				dist_threshold *= 0.9;
			}
		}

		return res;
	}

	bool validNewLocation(const typename TTree<DomCtr, TDim>::TPointD & location, double distance_threshold)
	{
		// generated point must be a certain distance away to all terminal points of the same tree
		for(unsigned int term_index : myIterator->myVectTerminals)
		{
			if((myIterator->myVectSegments[term_index].myCoordinate - location).norm() < distance_threshold)
			{
				return false;
			}
		}

		// generated point must be a certain distance away from every tree segments
		for(tree_vector_iterator it = myTrees.begin(); it < myTrees.end(); it++)
		{
			for(unsigned int i = 0; i < it->myVectSegments.size(); i++)
			{
				if(it->getProjDistance(i, location) < it->myVectSegments[i].myRadius)
				{
					return false;
				}
			}
		}

		// Other conditions for location validity can be placed here
		// maybe a template policy pattern would make this convenient

		return true;	// if we reach this, conditions are met
	}

	void replaceCurrentTree(const TTree<DomCtr, TDim> & tree)
	{
		*myIterator = tree;
	}

	bool isIntersectingTrees(const typename TTree<DomCtr, TDim>::TPointD & ptA,
							 const typename TTree<DomCtr, TDim>::TPointD & ptB,
							 double r, unsigned int seg_index)
	{
		//Get useful points
		unsigned int parent_ind = myIterator->myVectParent[seg_index];
		auto seg_start = myIterator->myVectSegments[parent_ind].myCoordinate;
		auto seg_end = myIterator->myVectSegments[seg_index].myCoordinate; 

		for(tree_vector_iterator it = myTrees.begin(); it != myTrees.end(); it++)
		{
			unsigned int LC_ind;	// variable reused for left child of considered segment
			unsigned int RC_ind;	// variable reused for right child of considered segment

			// if it's the same tree, the new segments are allowed to intersect their parent and children
			if(it == myIterator)
			{
				// middle segment
				LC_ind = it->myVectChildren[parent_ind].first;
				RC_ind = it->myVectChildren[parent_ind].second;
				if(it->isIntersectingTree(ptB, seg_start, r, std::make_tuple(LC_ind, RC_ind, parent_ind)))
				{
					return true;
				}

				// sibling segment (first child of middle segment)
				LC_ind = it->myVectChildren[seg_index].first;
				RC_ind = it->myVectChildren[seg_index].second;
				if(it->isIntersectingTree(ptB, seg_end, r, std::make_tuple(LC_ind, RC_ind, seg_index)))
				{
					return true;
				}

				// added segment (second child of middle segment)
				// since it's the added segment, it has no children
				LC_ind = -1;
				RC_ind = -1;
				if(it->isIntersectingTree(ptB, ptA, r, std::make_tuple(LC_ind, RC_ind, seg_index)))
				{
					return true;
				}
			}
			else
			{
				// other tree cases
				// no intersection is allower with any segment
				if(it->isIntersectingTree(ptB, seg_start, r, std::make_tuple(-1,-1,-1))
					|| it->isIntersectingTree(ptB, seg_end, r, std::make_tuple(-1,-1,-1))
				    || it->isIntersectingTree(ptB, ptA, r, std::make_tuple(-1,-1,-1)) )
				{
					return true;
				}
			}
		}

		return false;
	}

	void expansionSummary(bool verbose)
	{
		bool fully_expanded = true;
		for(std::size_t i = 0; i < myTrees.size(); i++)
		{
			fully_expanded = fully_expanded && (myTrees[i].my_NTerm == myTrees[i].myKTerm);
		}

		if(!fully_expanded)
		{
			DGtal::trace.warning() << "All seeds not found due to too large domain constraints :" << std::endl;

			for(std::size_t i = 0; i < myTrees.size(); i++)
			{
				DGtal::trace.warning() << "Tree " << i << " : " 
					<< myTrees[i].myKTerm << " / " << myTrees[i].my_NTerm << "." << std::endl;
			}
		}

		if(verbose)
		{
			for(std::size_t i = 0; i < myTrees.size(); i++)
			{
				std::cout << "Tree " << i << " ====> Aperf=" 
					<< myTrees[i].myRsupp * myTrees[i].myRsupp * myTrees[i].my_NTerm * M_PI
					<< " == " << myTrees[i].my_aPerf << std::endl;
			}
		}
	}

	unsigned int attemptsSum()
	{
		unsigned int res = 0;
		for(unsigned int attempts : myExpansionAttempts)
		{
			res += attempts;
		}

		return res;
	}

	unsigned int NTermsSum()
	{
		unsigned int res = 0;
		for(const TTree<DomCtr, TDim> & tree : myTrees)
		{
			res += tree.my_NTerm;
		}
		
		return res;
	}

private:
	tree_vector myTrees;
	tree_vector_iterator myIterator;
	std::vector<unsigned int> myExpansionAttempts;	// an attempt can be succesful or not; counts every try

	void initializeFirstSegments()
	{
		//
	}
};