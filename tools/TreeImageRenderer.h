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
	TreeImageRenderer(const unsigned int width,
					  const std::string & radii_filename,
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
void compBB(std::vector< typename TreeImageRenderer<TDim>::TPointD > &points,
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