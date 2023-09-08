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
#include <memory>

#include "DGtal/base/Common.h"

#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/ITKReader.h"
#include "DGtal/io/readers/VolReader.h"

#include "svg_elements.h"
#include "GeomHelpers.h"

// Aliases
template <int TDim>
using Image = DGtal::ImageContainerBySTLVector< Domain< Space< PointI<TDim> > >, double>;


///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////          TreeImageRenderer<TDim>              ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


template <int TDim>
class TreeImageRenderer
{
private:
	typedef Space< PointI<TDim> > TSpace;

public:
	// Nested classes

	struct Segment
	{
		// Distal point of the segment.
		unsigned int myDistalIndex;
		// Proximal point of the segment.
		unsigned int myProximalIndex;
		// hydrodynamc registance (R star)
		double myResistance = 0.00;
		// flow (Qi)
		double myFlow = 0.00;
	};


	struct ArteryTree
	{
		std::vector<Segment> mySegments;
		std::vector< PointD<TDim> > myPoints;
		std::vector<double> myRadii;			// Radii associated with the vertices
	};

	struct OrganDomain
	{
		OrganDomain() : myDomainMask(Domain<TSpace>()), isDefined(false)
		{
		}

		OrganDomain(const std::string & domain_filename);

		Image<TDim> myDomainMask;
		bool isDefined;
	};

	// Constructor

	/**
	 * @brief Constructor, the ArteryTree data vectors.
	 * @param radii_filename The name of the file containing radii data.
	 * @param vertices_filename The name of the file containing radii data.
	 * @param edges_filename The name of the file containing radii data.
	 **/
	TreeImageRenderer(const std::string & radii_filename,
					  const std::string & vertices_filename,
					  const std::string & edges_filename);


	/**
	 * @brief Constructor, the ArteryTree data vectors.
	 * @param radii_filename The name of the file containing radii data.
	 * @param vertices_filename The name of the file containing radii data.
	 * @param edges_filename The name of the file containing radii data.
	 * @param organ_dom_filename The image file of the organ domain.
	 **/
	TreeImageRenderer(const std::string & radii_filename,
					  const std::string & vertices_filename,
					  const std::string & edges_filename,
					  const std::string & organ_dom_filename);

	// Methods

	/**
	 * @brief Populates the artery tree vectors with the specified files.
	 * @param radii_filename The name of the file containing radii data.
	 * @param vertices_filename The name of the file containing radii data.
	 * @param edges_filename The name of the file containing radii data.
	 **/
	void importTreeData(const std::string & radii_filename,
						const std::string & vertices_filename,
						const std::string & edges_filename);


	/**
     * @brief Computes a render of myTree.
     * @brief The values of the pixels inside the segments depend on the flow value.
     * @param width The width in pixels of the output Image<TDim>.
     * @returns a Image<TDim> object.
     **/
	Image<TDim> flowRender(unsigned int width);


	Image<TDim> flowRender2(unsigned int width);

	/**
     * @brief Exports a 2D SVG animation of the construction of myTree.
     * @param filename The file in which the animation will be written.
     * @param duration The total duration of the animation
     **/
	void animationRender(const std::string & filename, int duration);


	/**
     * @brief Computes the skeleton of myTree.
     * @param width The width in pixels of the output Image<TDim>.
     * @returns a Image<TDim> object.
     **/
	Image<TDim> skeletonRender(unsigned int width);


	/**
     * @brief Computes the realistic render of myTree.
     * @param width The width in pixels of the output Image<TDim>.
     * @returns a Image<TDim> object.
     **/
	Image<TDim> realisticRender(double sigma, unsigned int width = 0);


private:
	/**
     * @brief Computes the image size given its desired width and the margins thickness.
     * @brief Will scale the data of myTree to fit the width, the height is controlled by the bounding box of myTree's points
     * @param width the total desired width in pixels, including margins, of the image.
     * @param margin_thickness the margin thickness in pixels.
     * @returns a Image<TDim> object of the desired width
     **/
	Image<TDim> createDomainImage(unsigned int width, unsigned int margin_thickness);

	// Member variables
	ArteryTree myTree;
	OrganDomain myOrganDomain;
};


///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////               Other functions                 ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////



/**
 * @brief Exports the 2D/3D image to the specified file.
 * @param image The image to export.
 * @param filename The file in which the image will be written.
 **/
template<int TDim>
void saveRender(const Image<TDim> & image,
				const std::string & filename);


/**
 * @brief Computes the upper and lower bounds of a vector of points.
 * @param[in] points The vector of points.
 * @param[out] upperbound The upperbound of the points.
 * @param[out] lowerbound The lowerbound of the points.
 **/
template<int TDim>
void compBB(const std::vector< PointD<TDim> > &points,
			PointD<TDim> &upperbound, 
			PointD<TDim> &lowerbound);


/**
 * @brief Projects ptC onto the line defined by ptA and ptB.
 * @param[in] ptA, ptB The points defining the segment.
 * @param[in] ptC The point to be projected.
 * @param[out] ptP The result of the projection.
 * @returns Whether the projection ptP belongs to the segment or not.
 **/
template<int TDim>
bool projectOnStraightLine(const PointD<TDim>& ptA,
						   const PointD<TDim> & ptB,
						   const PointD<TDim> & ptC,
						   PointD<TDim> & ptP);


/**
 * @brief Initializes a SVG::Line and its animations, if the pointer is nullptr (considered not already initialized)
 * @param proximal The starting point of the line
 * @param distal The ending point of the line
 * @param radius The line radius
 * @param color The color of the line
 * @param duration The duration of the animation
 * @param[out] line_ptr The line to initialize
 * @returns true if a SVG::Line was initialized.
 **/
bool initializeSVGLine(const PointD<2> & proximal,
					   const PointD<2> & distal, 
					   double radius,
					   const SVG::Color & color,
					   int duration,
					   std::shared_ptr<SVG::AnimatedElement> & line_ptr);


/**
 * @brief Draws a 1-pixel width line between 2 points, using Bresenham algorithm.
 * @brief It should always draw at least one pixel, even when p0 == p1
 * @param image The image where the line will be drawn (modified by the function).
 * @param p0 The first point of the line.
 * @param p1 The second point of the line.
 * @param value The value (understand color) of the line.
 **/
template<int TDim>
void drawBresenhamLine(Image<TDim> & image,
					   PointI<TDim> p0,
					   const PointI<TDim> & p1,
					   const double value)
{
	// check if line is within the image 
	if( !image.domain().isInside(p0) || !image.domain().isInside(p1) )
	{
		std::cout << "Points defining the line are not all in the image domain." << std::endl;
		return;
	}

	// Bresenham algorithm, adapted for n dimensions

	std::vector<int> deltas;				// contains dx, dy ...
	std::vector<int> increments;
	std::size_t driving_axis = 0;

	for(std::size_t i = 0; i < TDim; i++)
	{
		deltas.emplace_back(p1[i] - p0[i]);

		if(deltas[i] < 0)					// increment is 1 for positive slope, -1 otherwise
		{
			increments.emplace_back(-1);
			deltas[i] *= -1;				// from this point on, deltas values are positive
		}
		else
		{
			increments.emplace_back(1);
		}

		if(deltas[i] > deltas[driving_axis])
		{
			driving_axis = i;
		}
	}

	std::vector<int> diffs;					// errors compared to exact line
	for(std::size_t i = 0; i < TDim; i++)
	{
		diffs.emplace_back(2*deltas[i] - deltas[driving_axis]);
	}

	while(p0[driving_axis] - increments[driving_axis] != p1[driving_axis])
	{
		// draw point
		image.setValue(p0, value);
		
		// update current point and diffs
		for(std::size_t i = 0; i < TDim; i++)
		{
			// if i is the driving axis, we increment the coordinate no matter what
			if (i == driving_axis)
			{
				p0[i] += increments[i];
				continue;
			}

			// for the coordinates that are not the driving axis
			if(diffs[i] > 0)
			{
				p0[i] += increments[i];
				diffs[i] += 2 * (deltas[i] - deltas[driving_axis]);
			}
			else
			{
				diffs[i] += 2 * deltas[i];	
			}
		}
	}
}


/**
 * @brief Brings all values of the Image into a range (default range is [0.0, 1.0].
 * @param image The image to normalize.
 * @param lb The lower bound of the range.
 * @param ub The upper bound of the range.
 **/
template <int TDim>
void normalizeImageValues(Image<TDim> & image, double lb = 0.0, double ub = 1.0)
{
	auto min_val = std::min_element(image.constRange().begin(), image.constRange().end());
	auto max_val = std::max_element(image.constRange().begin(), image.constRange().end());

	if(*min_val == *max_val || lb == ub)
	{
		for(auto it = image.begin(); it != image.end(); it++)
		{
			*it = lb;
		}
	}
	else
	{
		double k = (ub - lb) / (*max_val - *min_val); 		// != 0 and != inf

		for(auto it = image.begin(); it != image.end(); it++)
		{
			*it = lb + (*it - *min_val) * k;
		}
	}
	
	return;
}


/**
 * @brief Blends two images together.
 * @brief Depending on the implementation, the order of the two input images may matter.
 * @brief All three images (input and output) MUST have the same domain.
 * @param[in] img1 The first image.
 * @param[in] img2 The second image.
 * @param[in] blendmode The blending logic for each resulting pixel. MUST have signature (double)(double, double).
 * @param[out] imgres The resulting image.
 **/
template <int TDim>
void imageBlend(const Image<TDim> & img1,
				const Image<TDim> & img2,
				double (*blendmode)(double, double),
				Image<TDim> & imgres)
{
	// check if all images have the same dimensions
	if(img1.domain().lowerBound() != img2.domain().lowerBound()
	   || img1.domain().lowerBound() != imgres.domain().lowerBound()
	   || img1.domain().upperBound() != img2.domain().upperBound()
	   || img1.domain().upperBound() != imgres.domain().upperBound())
	{
		return;
	}

	for(const PointI<TDim> & p : img1.domain())
	{
		imgres.setValue(p, blendmode(img1(p), img2(p)));
	}
}