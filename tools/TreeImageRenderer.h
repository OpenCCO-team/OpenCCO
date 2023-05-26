#pragma once

#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include "GeomHelpers.h"

typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, double> TImage2D;
typedef DGtal::Z2i::Point TPoint;
typedef DGtal::PointVector<2, double> TPointD;
typedef DGtal::Z2i::Domain TDomain;

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

class TreeImageRenderer
{
public:
	// Constructor
	TreeImageRenderer(unsigned int width);

	// Methods
	void importTreeData();

	/**
     * Computes the image size given its desired width and the margins thickness
     * @param width the total desired width in pixels, including margins, of the image
     * @param margin_thickness the margin thickness in pixels
     **/
	void setImageSize(unsigned int width, unsigned int margin_thickness);

	void createDistanceAlphaChannel();

private:
	// Member variables
	TImage2D myBackground;
	TImage2D myTreeImage;
	TImage2D myDistanceMap;
	ArteryTree myTree;
	TDomain myDomain;
};

// Other useful fonctions

/**
 * Computes the upper and lower bounds of a vector of points
 * @param points the vector of points
 * @param upperbound the upperbound of the points (modified by this function)
 * @param lowerbound the lowerbound of the points (modified by this function)
 **/
void compBB(std::vector<TPointD> &points, TPointD &upperbound, TPointD &lowerbound);



/**
 * Projects ptC onto the line defined by ptA and ptB
 * @param ptA 1st point of the segment
 * @param ptB 2nd point of the segment
 * @param ptC the point to be projected
 * @param ptP the result of the projection (modified by this function)
 * @return true if the projection belongs to the segment
 **/
bool projectOnStraightLine(const TPointD & ptA,
						   const TPointD & ptB,
						   const TPointD & ptC,
						   TPointD & ptP);



/**
 * Creates a subdomain from a domain and desired lower and upper bounds.
 * Upper and lower bounds may be outside the domain, in other words the 
 * @param dom1 1st domain
 * @param dom2 2nd domain
 * @return the domain resulting from the intersection
 **/
TDomain domainIntersect(TDomain &dom1, TDomain &dom2);


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