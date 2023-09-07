#pragma once

#include <iostream>

#include "XmlHelpers.h"
#include "CoronaryArteryTree.h"


/**
 * The goal of CohabitingTrees template class is to manage multiple trees in the same domain.
 * It can check the state of each tree, export them, etc
 * 
 * Will manage which tree needs to be expanded next.
 */
template<class DomCtr, int TDim>
class CohabitingTrees
{
	typedef std::vector< CoronaryArteryTree<DomCtr, TDim> > tree_vector;
	typedef typename std::vector< CoronaryArteryTree<DomCtr, TDim> >::iterator tree_vector_iterator;

public:

	/**
     * @brief Constructor.
     * @brief It initializes the desired number of trees and their first segments.
     * @param aPerf Surface of the perfusion (common to all trees).
     * @param nTerm_vec Number of terminal segments, governs the number of trees
     * 		(each tree has a specified number of terminal segments).
     * @param dom_ctr The domain controller, common for each tree.
     * @param radius The radius of the root segment.
     **/
	CohabitingTrees(double aPerf, const std::vector<unsigned int> & nTerm_vec, double gamma,
					std::vector< PointI<TDim> > starting_points,
					DomCtr & dom_ctr, double radius = 1.0)
	{
		int n = nTerm_vec.size();

		for(std::size_t i = 0; i < n; i++)
		{
			myTrees.emplace_back(aPerf, nTerm_vec[i], dom_ctr, radius);
			myExpansionAttempts.push_back(1);	// inital value is one because the first segment counts as a succesful attempt
		}

		if(starting_points.empty())
		{
			starting_points = firstN_CandidatePoints<DomCtr, TDim>(dom_ctr, n);
		} 


		// initialize first segment of trees so that they're not overlapping
		for(std::size_t i = 0; i < n; i++)
		{
			myTrees[i].myVectSegments[0].myCoordinate = starting_points[i];
			myTrees[i].myVectSegments[1].myCoordinate = starting_points[i] - 0.25 * starting_points[i];
		}

		myIterator = myTrees.begin();
	}

	// Methods

	/**
	 * @brief Method to set myIterator to the next tree.
	 * @brief Right now it increments myIterator, and if it reaches the end makes it loop back.
	 * @brief An upgrade could be to use a template policy to have maximum flexibility.
	 **/
	void nextTree()
	{
		myIterator++;

		// check if myIterator is valid
		if(myIterator < myTrees.begin() || myIterator >= myTrees.end())
		{
			myIterator = myTrees.begin();
		}
	}


	/**
	 * Creates a copy of the tree currently pointed to by myIterator.
	 * @returns a CoronaryArteryTree<DomCtr, TDim>.
	 **/
	CoronaryArteryTree<DomCtr, TDim> getCurrentTreeCopy()
	{
		return *myIterator;
	}


	/**
	 * @brief Shortcut method to get the radius from the segment index of the current tree.
	 * @returns The radius of the segment.
	 */
	double getSegmentRadius(unsigned int seg_index)
	{
		return myIterator->myVectSegments[seg_index].myRadius;
	}


	/**
	 * @brief Finds the n nearest segment from the current tree closest to a point in the domain.
	 * @param p The coordinates of the point.
	 * @returns A vector containing the indices of the closest segments.
	 * 
	 */
	std::vector<unsigned int> getNeighbors(const PointD<TDim> & p)
	{
		return myIterator->getN_NearestSegments(p, myIterator->myNumNeighbor);
	}


	/**
	 * @brief Increments the tally of attemps made to expand the current tree.
	 */
	void incrementAttempt()
	{
		std::size_t i = myIterator - myTrees.begin();
		myExpansionAttempts[i]++;
	}


	/**
	 * @brief Evaluates if the trees have been expanded as much as possible given constraints.
	 * @returns whether or not the expansion is finished 
	 **/
	bool expansionFinished()
	{
		bool res = true;

		for(std::size_t i = 0; i < myTrees.size(); i++)
		{
			res = res && (myExpansionAttempts[i] >= myTrees[i].my_NTerm);
		}

		return res;
	}


	/**
	 * @brief Randomly generates a random point within the domain controller.
	 * @brief The point respects a distance constraint.
	 * @param nb_trials The number of point to be generated before relaxing the distance constraint.
	 * @returns a point.
	 **/
	PointD<TDim> generateNewLocation(unsigned int nb_trials)
	{
		PointD<TDim> res;
		double dist_threshold = myIterator->getDistanceThreshold();
		
		unsigned int i = 0;
		bool valid = false;
		while(!valid)
		{
			res = myIterator->myDomainController().randomPoint();

			valid = validNewLocation(res, dist_threshold);

			// if we reach end of the loop before finding the new location,
			// loop again with less restrictive distance threshold
			i++;			
			if(i == nb_trials)
			{
				i = 0;
				dist_threshold *= 0.9;
			}
		}

		return res;
	}


	/**
	 * @brief Method that verifies that a point is valid considering a distance constraint.
	 * @param location The point to evaluate.
	 * @param distance_threshold The distance constraint.
	 * @returns Whether or not the point is valid
	 **/
	bool validNewLocation(const PointD<TDim> & location, double distance_threshold)
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


	/**
	 * @brief replaces the current tree (pointed to by myIterator) with another.
	 * @param tree The new tree.
	 **/
	void replaceCurrentTree(const CoronaryArteryTree<DomCtr, TDim> & tree)
	{
		*myIterator = tree;
	}

	/**
	 * @brief Evaluates if a segment is intersecting with other segments of any trees.
	 * @param p_added The end point of the segment.
	 * @param p_bifurcation The starting point of the segment (common with it's parent and sibling).
	 * @param r The radius of the segment.
	 * @param seg_index The index of the parent segment.
	 * @returns Whether the segment is intersecting with any other segment.
	 **/
	bool isIntersectingTrees(const PointD<TDim> & p_added,
							 const PointD<TDim> & p_bifurcation,
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
				if(it->isIntersectingTree(p_bifurcation, seg_start, r, std::make_tuple(LC_ind, RC_ind, parent_ind)))
				{
					return true;
				}

				// sibling segment (first child of middle segment)
				LC_ind = it->myVectChildren[seg_index].first;
				RC_ind = it->myVectChildren[seg_index].second;
				if(it->isIntersectingTree(p_bifurcation, seg_end, r, std::make_tuple(LC_ind, RC_ind, seg_index)))
				{
					return true;
				}

				// added segment (second child of middle segment)
				// since it's the added segment, it has no children
				LC_ind = -1;
				RC_ind = -1;
				if(it->isIntersectingTree(p_bifurcation, p_added, r, std::make_tuple(LC_ind, RC_ind, seg_index)))
				{
					return true;
				}
			}
			else
			{
				// other tree cases
				// no intersection is allowed with any segment
				if(it->isIntersectingTree(p_bifurcation, seg_start, 8*r, std::make_tuple(-1,-1,-1))
					|| it->isIntersectingTree(p_bifurcation, seg_end, 8*r, std::make_tuple(-1,-1,-1))
					|| it->isIntersectingTree(p_bifurcation, p_added, 8*r, std::make_tuple(-1,-1,-1)) )
				{
					return true;
				}
			}
		}

		return false;
	}


	/**
	 * @brief Displays in the console a summary of the expansion.
	 * @brief If the expansion was not complete (in the sense of every tree having the desired amount of terminal segments),
	 * @brief displays the details of each tree.
	 * @param verbose To display or not more data.
	 **/
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


	/**
	 * @brief Sums the attempts of expansion made for each tree.
	 * @returns the sum of the attempts.
	 **/
	unsigned int attemptsSum() const
	{
		unsigned int res = 0;
		for(unsigned int attempts : myExpansionAttempts)
		{
			res += attempts;
		}

		return res;
	}


	/**
	 * @brief Sums the desired number of terminal segments for each tree.
	 * @returns the total.
	 **/
	unsigned int NTermsSum() const
	{
		unsigned int res = 0;
		for(const CoronaryArteryTree<DomCtr, TDim> & tree : myTrees)
		{
			res += tree.my_NTerm;
		}
		
		return res;
	}


	/**
	 * @brief If the dimension is 2, creates an svg file representing the trees in the domain.
	 **/
	void exportTreesDisplays()
	{
		if(TDim == 2)	// Can only create a 2D display for 2D generated trees.
		{
			// DGtal board to draw the different elements
			DGtal::Board2D board;
			board.setLineCap(LibBoard::Shape::LineCap::RoundCap);
			board.setUnit(0.1, LibBoard::Board::UCentimeter);

			// Colormaps, one for each of the tree so that they're easy to recognize of the svg
			std::vector< DGtal::GradientColorMap<float> > colormaps;
			for(std::size_t i = 0; i < myTrees.size(); i++)
			{
				colormaps.emplace_back(std::log(myTrees[i].my_qTerm), std::log(myTrees[i].my_qPerf));

				// each tree has a different color, we use the hue
				double h = i * 360.0 / myTrees.size();
				DGtal::Color c_low;
				DGtal::Color c_high;
				c_low.setFromHSV(h, 0.10, 0.90);		// desaturated and lighter
				c_high.setFromHSV(h, 0.95, 0.784);		// saturated and a bit darker

				colormaps[i].addColor(c_low);
				colormaps[i].addColor(c_high);
			}

			// Draw the trees
			for(std::size_t i = 0; i < myTrees.size(); i++)
			{
				for(const typename CoronaryArteryTree<DomCtr, TDim>::Segment & s : myTrees[i].myVectSegments)
				{
					if(s.myIndex == 0)
						continue;

					board.setPenColor(colormaps[i](std::log(s.myFlow)));

					PointD<TDim> proximal = myTrees[i].myVectSegments[myTrees[i].myVectParent[s.myIndex]].myCoordinate;
					PointD<TDim> distal = s.myCoordinate;

					board.setLineWidth(s.myRadius);
					board.drawLine(distal[0], distal[1], proximal[0], proximal[1]);
				}
			}

			board.saveSVG("cohabiting_trees.svg");
		}
	}


	/**
	 * @brief Creates an xml file for each tree, like the classic generation would.
	 * @brief Default exports are like "tree_<TDim>_<n>.xml" where n is the index of the tree.
	 * @brief @param filename the prefix of the xml files : <filename>_<n>.xml
	 **/
	void writeTreesToXML(std::string filename)
	{
		for(std::size_t i = 0; i < myTrees.size(); i++)
		{
			std::string xml_filename = filename + "_" + std::to_string(i) + ".xml";
			XMLHelpers::writeTreeToXml(myTrees[i], xml_filename.c_str());
		}
	}


private:
	// Member variables
	tree_vector myTrees;							// Tree list
	tree_vector_iterator myIterator;				// current tree
	std::vector<unsigned int> myExpansionAttempts;	// An attempt can be succesful or not; counts every try
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
std::vector< PointI<3> > evenlySpreadPoints3D(const DomCtr & domain, unsigned int n)
{
	std::cout << domain.lowerBound()[0] << " " << domain.lowerBound()[1] << " " << domain.lowerBound()[2] << std::endl;
	std::cout << domain.upperBound()[0] << " " << domain.upperBound()[1] << " " << domain.upperBound()[2] << std::endl;
	std::cout << domain.myCenter[0] << " " << domain.myCenter[1] << " " << domain.myCenter[2] << std::endl << std::endl;

	std::vector< PointD<3> > base_points;
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
		for(const PointD<3> & p : base_points)
		{
			within_r1 = within_r1 && domain.isInside(PointI<3>(domain.myCenter + r1 * p));
			within_r2 = within_r2 && domain.isInside(PointI<3>(domain.myCenter + r2 * p));
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
		for(const PointD<3> & p : base_points)
		{
			within_r_mean = within_r_mean && domain.isInside(PointI<3>(domain.myCenter + r_mean * p));
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

	std::vector< PointI<3> > res;

	for(const PointD<3> & p : base_points)
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
std::vector< PointI<2> > evenlySpreadPoints2D(const DomCtr & domain, unsigned int n)
{
	std::vector< PointD<2> > base_points;
	const double delta_theta = 2 * M_PI / n;	// var to avoid doing operation each loop

	for(std::size_t i = 0; i < n; i++)
	{
		// angle
		double theta = i * delta_theta;
		base_points.emplace_back(cos(theta), sin(theta));
	}

	return std::vector<PointI<2>>();
}