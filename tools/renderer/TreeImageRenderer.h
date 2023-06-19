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

#include "svg_elements.h"


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
	};


	struct ArteryTree
	{
		std::vector<Segment> mySegments;
		std::vector<TPointD> myPoints;
		std::vector<double> myRadii;			// Radii associated with the vertices
	};

	// Constructor

	/**
	 * @brief Constructor, initializes the size of the image and the artery tree data vectors.
	 * @param width The desired width of the image.
	 * @param radii_filename The name of the file containing radii data.
	 * @param vertices_filename The name of the file containing radii data.
	 * @param edges_filename The name of the file containing radii data.
	 **/
	TreeImageRenderer(const std::string & radii_filename,
					  const std::string & vertices_filename,
					  const std::string & edges_filename);

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
     * @brief Computes the image size given its desired width and the margins thickness.
     * @param width the total desired width in pixels, including margins, of the image.
     * @param margin_thickness the margin thickness in pixels.
     **/
	void setImageSize(unsigned int width, unsigned int margin_thickness);


	/**
     * @brief Computes the distance tranform of the artery tree.
     * @brief The result is stored in myDistanceMap.
     **/
	void createDistanceMap();
	

	/**
     * @brief Computes the image of the artery tree.
     * @brief The result is stored in myTreeImage.
     **/
	void createTreeImage();


	/**
     * @brief Exports the 2D/3D image to the specified file.
     * @param filename The file in which the image will be written.
     **/
	void saveRender(const std::string & filename);

	/**
     * @brief Exports a 2D SVG animation of the construction of the ArteryTree.
     * @param filename The file in which the animation will be written.
     * @param duration The total duration of the animation
     **/
	void treeConstructionAnimation(const std::string & filename, int duration);

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
 * @brief Computes the upper and lower bounds of a vector of points.
 * @param[in] points The vector of points.
 * @param[out] upperbound The upperbound of the points.
 * @param[out] lowerbound The lowerbound of the points.
 **/
template<int TDim>
void compBB(const std::vector< typename TreeImageRenderer<TDim>::TPointD > &points,
			typename TreeImageRenderer<TDim>::TPointD &upperbound, 
			typename TreeImageRenderer<TDim>::TPointD &lowerbound);



/**
 * @brief Projects ptC onto the line defined by ptA and ptB.
 * @param[in] ptA, ptB The points defining the segment.
 * @param[in] ptC The point to be projected.
 * @param[out] ptP The result of the projection.
 * @returns Whether the projection ptP belongs to the segment or not.
 **/
template<int TDim>
bool projectOnStraightLine(const typename TreeImageRenderer<TDim>::TPointD & ptA,
						   const typename TreeImageRenderer<TDim>::TPointD & ptB,
						   const typename TreeImageRenderer<TDim>::TPointD & ptC,
						   typename TreeImageRenderer<TDim>::TPointD & ptP);


/**
 * @brief Initializes a SVG::Line and its animations, if the pointer is nullptr (considered not already initialized)
 * @param proxital The starting point of the line
 * @param distal The ending point of the line
 * @param radius The line radius
 * @param color The color of the line
 * @param duration The duration of the animation
 * @param[out] line_ptr The line to initialize
 * @returns true if a SVG::Line was initialized.
 **/
bool initializeSVGLine(const TreeImageRenderer<2>::TPointD & proxital,
					   const TreeImageRenderer<2>::TPointD & distal, 
					   double radius,
					   const SVG::Color & color,
					   int duration,
					   std::shared_ptr<SVG::AnimatedElement> & line_ptr);


/**
 * @brief Draws a 1-pixel width line between 2 points, using Bresenham algorithm.
 * @param image The image where the line will be drawn (modified by the function).
 * @param p0 The first point of the line.
 * @param p1 The second point of the line.
 **/
template<int TDim>
void drawBresenhamLine(typename TreeImageRenderer<TDim>::TImage & image,
					   typename TreeImageRenderer<TDim>::TPoint p0,
					   const typename TreeImageRenderer<TDim>::TPoint & p1)
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
		image.setValue(p0, 255);
		
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