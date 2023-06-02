#pragma once

#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <limits>
#include <sstream>

#include "DGtal/base/Common.h"

#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/helpers/StdDefs.h"


///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                   Segment                     ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


struct Segment
{
	// Distal point of the segment.
	unsigned int myDistalIndex;
	// Proxital point of the segment.
	unsigned int myProxitalIndex;
	// hydrodynamc registance (R star)
	double myResistance = 0.00;
	// flow (Qi)
	double myFlow = 0.00;

	// Operator overloading
	// friend keyword has no effect since all member variable are public
	// but I've heard it's good practice
	friend bool operator<(const Segment & s1, const Segment & s2);
};


///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////            TreeImageRenderer                  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


template <int TDim>
class TreeImageRenderer
{
public:
	typedef DGtal::SpaceND<TDim, int>   SpaceCT;
	typedef DGtal::PointVector<TDim, int> TPoint;
	typedef DGtal::PointVector<TDim, double> TPointD;
	typedef DGtal::HyperRectDomain<SpaceCT> TDomain;
	typedef DGtal::ImageContainerBySTLVector<TDomain, double> TImage;

	// Nested classes

	struct ArteryTree
	{
		std::vector<Segment> mySegments;
		std::vector<TPointD> myPoints;
		std::vector<double> myRadii;			// Radii associated with the vertices

		// Method


		/**
	     * Sorts mySegments, myPoints and myRaddi based on the positions of the points.
	     * The points are compared on their first coordinate minus their radius
	     **/
		void sort();
	};

	// Constructor
	TreeImageRenderer(const unsigned int width,
					  const std::string & radii_filename,
					  const std::string & vertices_filename,
					  const std::string & edges_filename);

	// Methods
	void importTreeData(const std::string & radii_filename,
						const std::string & vertices_filename,
						const std::string & edges_filename);

	/**
     * Computes the image size given its desired width and the margins thickness
     * @param width the total desired width in pixels, including margins, of the image
     * @param margin_thickness the margin thickness in pixels
     **/
	void setImageSize(unsigned int width, unsigned int margin_thickness);

	void createDistanceMap();
	
	void createTreeImage();

	void saveRender(const std::string & filename);

	const TImage & distanceMap() const;

	const TImage & treeImage() const;

	static const int myDim = TDim;

private:
	// Member variables
	TImage myBackground;
	TImage myTreeImage;
	TImage myDistanceMap;
	ArteryTree myTree;
	TDomain myDomain;
};


///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////               Other functions                 ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////



/**
 * Computes the upper and lower bounds of a vector of points
 * @param points the vector of points
 * @param upperbound the upperbound of the points (modified by this function)
 * @param lowerbound the lowerbound of the points (modified by this function)
 **/
template<int TDim>
void compBB(std::vector< typename TreeImageRenderer<TDim>::TPointD > &points,
			typename TreeImageRenderer<TDim>::TPointD &upperbound, 
			typename TreeImageRenderer<TDim>::TPointD &lowerbound);



/**
 * Projects ptC onto the line defined by ptA and ptB
 * @param ptA 1st point of the segment
 * @param ptB 2nd point of the segment
 * @param ptC the point to be projected
 * @param ptP the result of the projection (modified by this function)
 * @return true if the projection belongs to the segment
 **/
template<int TDim>
bool projectOnStraightLine(const typename TreeImageRenderer<TDim>::TPointD & ptA,
						   const typename TreeImageRenderer<TDim>::TPointD & ptB,
						   const typename TreeImageRenderer<TDim>::TPointD & ptC,
						   typename TreeImageRenderer<TDim>::TPointD & ptP);



// from https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(const std::vector<T> & vec,
										  Compare compare)
{
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(),
        [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
    return p;
}

// from https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
template <typename T>
std::vector<T> apply_permutation(const std::vector<T> & vec,
    							 const std::vector<std::size_t> & p)
{
    std::vector<T> sorted_vec(vec.size());
    std::transform(p.begin(), p.end(), sorted_vec.begin(),
        [&](std::size_t i){ return vec[i]; });
    return sorted_vec;
}