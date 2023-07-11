#pragma once

#include <iostream>
#include "CoronaryArteryTree.h"


template <class DomCtr, int TDim>
using TTree = CoronaryArteryTree<DomCtr, TDim>;

template <int TDim>
using TPoint = DGtal::PointVector<TDim, int>;

template <int TDim>
using TPointD = DGtal::PointVector<TDim, double>;

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

	double getSegmentRadius(unsigned int seg_index)
	{
		return myIterator->myVectSegments[seg_index].myRadius;
	}

	std::vector<unsigned int> getNeighbors(const TPointD<TDim> & p)
	{
		return myIterator->getN_NearestSegments(p, myIterator->myNumNeighbor);
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

	TPointD<TDim> generateNewLocation(unsigned int nb_trials)
	{
		TPointD<TDim> res;
		double dist_threshold = myIterator->getDistanceThreshold();
		
		unsigned int i = 0;
		while(i < nb_trials)
		{
			res = myIterator->myDomainController().randomPoint();

			if(validNewLocation(res, dist_threshold))
			{
				return res;
			}

			// if we reach end of the loop before finding the new location,
			// loop again with less restrictive distance threshold
			i++;			
			if(i == nb_trials)
			{
				i = 0;
				dist_threshold *= 0.9;
			}
		}
	}

	bool validNewLocation(const TPointD<TDim> & location, double distance_threshold)
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

	bool isIntersectingTrees(const TPointD<TDim> & ptA,
							 const TPointD<TDim> & ptB,
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

	void expansionSummary(bool verbose) const
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

	unsigned int attemptsSum() const
	{
		unsigned int res = 0;
		for(unsigned int attempts : myExpansionAttempts)
		{
			res += attempts;
		}

		return res;
	}

	unsigned int NTermsSum() const
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

	bool initializeFirstSegments(const std::vector< TPointD<TDim> > & starting_points)
	{
		// error if number of points != number of trees
		if(starting_points.size() != myTrees.size())
		{
			std::cout << "Amount of starting points is incorrect : " 
				<< starting_points.size() << " but should be " << myTrees.size() << std::endl;
			return false;
		}

		for(std::size_t i = 0; i < myTrees.size(); i++)
		{
			myTrees[i].myTreeCenter = myTrees[i].myDomainController().myCenter;
			myTrees[i].myVectSegments[0].myCoordinate = starting_points[i];
		}

		return true;
	}
};



/*
 * @brief returns a list of points of a 3D sphere.
 * @brief The points are placed with the golden spiral method to ensure they're evenly spread reliably.
 * @brief the sphere is centered on the DomainController center, and as big as possible.
 * @param domain The domain to fit the sphere.
 * @param n The amount of points to spread.
 * @returns a std::vector of points.
 */
template<class DomCtr>
std::vector< TPoint<3> > evenlySpreadPoints3D(const DomCtr & domain, unsigned int n)
{
	std::cout << domain.lowerBound()[0] << " " << domain.lowerBound()[1] << " " << domain.lowerBound()[2] << std::endl;
	std::cout << domain.upperBound()[0] << " " << domain.upperBound()[1] << " " << domain.upperBound()[2] << std::endl;
	std::cout << domain.myCenter[0] << " " << domain.myCenter[1] << " " << domain.myCenter[2] << std::endl << std::endl;

	std::vector< TPointD<3> > base_points;
	double epsilon = 0.36;				// offset for the indices, proven best for most n values : https://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/#more-3069
	double double_GR = (1 + sqrt(5));	// golden ratio times 2 (to save a multiplication by 2 later)

	for(std::size_t i = 0; i < n; i++)
	{
		// angles
		double phi = acos(1 - 2 * (i + epsilon) / n);
		double theta = M_PI * double_GR * i;

		base_points.emplace_back(cos(theta) * sin(phi),
						 sin(theta) * sin(phi),
						 cos(phi) );
	}

	// find the largest radius for which every point is within the domain
	// using a dichotomy
	double r1 = domain.myRadius;
	double r2 = 2 * r1;
	double precision = 0.5;

	// first we find r1 and r2 such that r2 = 2*r1 and all points with r1 are within the domain, and at least one point with r2 is not
	bool r_bounds_set = false;
	int DEBUG_I = 0;
	while(!r_bounds_set)
	{
		DEBUG_I++;
		bool within_r1 = true;			// whether all points with r1 are inside the domain
		bool within_r2 = true;			// whether all points with r2 are inside the domain
		for(const TPointD<3> & p : base_points)
		{
			within_r1 = within_r1 && domain.isInside(TPoint<3>(domain.myCenter + r1 * p));
			within_r2 = within_r2 && domain.isInside(TPoint<3>(domain.myCenter + r2 * p));
		}

		r_bounds_set = within_r1 && !within_r2;

		// adjust bounds if they don't satisfy criteria
		if(within_r1)
		{	
			if(within_r2)
			{
				r1 = r2;
				r2 *= 2;
			}
		}
		else
		{
			if(within_r2)
			{
				// r1 > r2 shouldn't ever happen, but managing this case anyway
				std::swap(r1, r2);
			}
			else
			{
				r2 = r1;
				r1 /= 2;
			}
		}
	}

	std::cout << r1 << " " << r2 << " found in " << DEBUG_I << " loops." << std::endl << std::endl;

	// dichotomy main loop
	DEBUG_I = 0;
	while(fabs(r1 - r2) > precision)
	{
		DEBUG_I++;
		double r_mean = (r1 + r2) / 2.0;

		double within_r_mean = true;
		for(const TPointD<3> & p : base_points)
		{
			within_r_mean = within_r_mean && domain.isInside(TPoint<3>(domain.myCenter + r_mean * p));
		}

		// modify r1 or r2, while keeping their property
		if(within_r_mean)
		{
			r1 = r_mean;
		}
		else
		{
			r2 = r_mean;
		}
	}

	std::cout << r1 << " " << r2 << " found in " << DEBUG_I << " loops." << std::endl;

	std::vector< TPoint<3> > res;

	for(const TPointD<3> & p : base_points)
	{
		res.emplace_back(p);
	}

	return res;
}



/*
 * @brief returns a list of points of a 2D sphere (i.e a circle).
 * @brief The points are uniformly spread around the edge of the circle.
 * @brief The circle is centered on the DomainController center, and as big as possible.
 * @param domain The domain to fit the sphere.
 * @param n The amount of points to spread.
 * @returns a std::vector of points.
 */
template<class DomCtr>
std::vector< TPoint<2> > evenlySpreadPoints2D(const DomCtr & domain, unsigned int n)
{
	std::vector< TPointD<2> > base_points;
	const double delta_theta = 2 * M_PI / n;	// var to avoid doing operation each loop

	for(std::size_t i = 0; i < n; i++)
	{
		// angle
		double theta = i * delta_theta;
		base_points.emplace_back(cos(theta), sin(theta));
	}

	return std::vector<TPoint<2>>();
}