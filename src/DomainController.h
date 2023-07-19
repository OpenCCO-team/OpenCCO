/**
 *  OpenCCO implementation 
 *  Copyright (C) 2023 B. Kerautret;  Phuc Ngo, N. Passat H. Talbot and C. Jaquet
 *
 *  This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 **/

#pragma once

#if defined(DOMAIN_CONTROLLER_RECURSE)
#error Recursive header files inclusion detected in DomainController.h
#else // defined(DOMAIN_CONTROLLER_RECURSE)
/** Prevents recursive inclusion of headers. */
#define DOMAIN_CONTROLLER_RECURSE

#if !defined DOMAIN_CONTROLLER_H
/** Prevents repeated inclusion of headers. */
#define DOMAIN_CONTROLLER_H


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/io/readers/GenericReader.h"
#include "GeomHelpers.h"


/**
 * Classes to control the construction domain of the coronary aretery tree.
 *  It is used during the reconstruction using the following points:
 *  - myCenter: used to orient the construction according the domain.
 *  - bool isInside(const TPoint &p):
 *    Test if a point is inside the domain or not (a point can be outside after the optmisation).
 *  - TPoint  randomPoint():
 *    To get a point inside the domaine used in the tree construction.
 *  - bool checkNoIntersectDomain(const TPoint &pt1, const TPoint &pt2):
 Usefull to ensure is a whole segment is inside the domain.
 *  - maxDistantPointFromBorder() const:
 *    usefull to determine a starting point and help to construct the tree.
 *  - TPoint firstCandidatePoint() const :
 *    use to initiate the reconstruction of the tree.
 **/


/**
 *  Domain controller based on a ball
 */
template<int TDim>
class CircularDomainCtrl
{
public:
	typedef enum {NO_UPDATE, UPDATED} UPDATE_RAD_TYPE ;
	UPDATE_RAD_TYPE myUpdateType = UPDATED;
	double myRadius {1.0};
	PointD<TDim> myCenter;
	
	// Constructor for ImplicitCirc type
	CircularDomainCtrl(){};
	
	CircularDomainCtrl(double radius, const PointD<TDim> & center)
	{
		myRadius = radius;
		myCenter = center;
	};
	
	virtual bool isInside(const PointD<TDim> & p) const
	{
		return (myCenter-p).norm() < myRadius;
	}
	
	PointD<TDim>
	randomPoint()
	{
		bool found = false;
		PointD<TDim> p;
		while(!found)
		{
			double ss = 0.0;
			for(unsigned int i = 0;i<TDim; i++ )
			{
				p[i] = ((double)rand() / RAND_MAX)*2.0*myRadius - myRadius;
			}

			p += myCenter;
			found = isInside(p);
		}
		return p;
	}

	/**
	 * Get the supported domain of the tree. By default it is defined from the circle center.
	 * If the domain if defined from a mask image, the center if computed from the imate center.
	 */
	PointD<TDim> getDomainCenter() const
	{
		return myCenter;
	}
	
	/**
	 * Check if the segment defined by two points intersect the domain.
	 *
	 * @param pt1 first point of the segment
	 * @param pt2  second point of the segment
	 */
	bool
	checkNoIntersectDomain(const PointD<TDim> & pt1, const PointD<TDim> & pt2)
	{
		return isInside(pt1) && isInside(pt2);
	}
	
	PointD<TDim>
	maxDistantPointFromBorder() const
	{
		return myCenter;
	}
	
	PointD<TDim>
	firstCandidatePoint() const
	{
		PointD<TDim> res;
		if (TDim == 2)
		{
			res[0] = 0.0;
			res[1] = myRadius;
		}
		else  if (TDim == 3)
		{
			res[0] = 0.0;
			res[1] = myRadius;
			res[2] = 0.0;
		}
		
		return res;
	}

	std::vector<std::vector< PointD<TDim> > >
	contours()
	{
		std::vector<std::vector< PointD<TDim> > > res;
		return res;
	}

	PointD<TDim>
	lowerBound()
	{
		PointD<TDim> p = PointD<TDim>::diagonal(myRadius*0.01);
		return myCenter - p;
	}

	PointD<TDim>
	upperBound()
	{
		PointD<TDim> p = PointD<TDim>::diagonal(myRadius*0.01);
		return myCenter + p;
	} 
};


template<int TDim>
class SquareDomainCtrl: public CircularDomainCtrl<TDim>{
	
public:
	// Constructor for ImplicitCirc type
	SquareDomainCtrl(){};
	
	SquareDomainCtrl(double radius,
					 const PointD<TDim> & center)
	{
		CircularDomainCtrl<TDim>::myRadius = radius;
		CircularDomainCtrl<TDim>::myCenter = center;
	};

	bool isInside(const PointD<TDim> & p) const
	{
		bool res = true;
		for (unsigned int i=0; i<TDim; i++)
		{
			res = res && (CircularDomainCtrl<TDim>::myCenter-p)[i] < CircularDomainCtrl<TDim>::myRadius;
		}
		return res;
	}
};


/**
 *  Domain controller based on a Image Mask
 */
template<int TDim>
class ImageMaskDomainCtrl
{
public:
	typedef DGtal::HyperRectDomain< Space< PointI<TDim> > > DomCT;
	typedef DGtal::ImageContainerBySTLVector< DomCT, unsigned char> Image;
	typedef DGtal::ImageContainerBySTLVector< DomCT, double> ImageD;
	typedef typename DGtal::DigitalSetSelector<DomCT, DGtal::BIG_DS + DGtal::HIGH_BEL_DS>::Type TDGset;
	typedef enum {NO_UPDATE, UPDATED} UPDATE_RAD_TYPE ;
	
	UPDATE_RAD_TYPE myUpdateType = NO_UPDATE;
	
	PointD<TDim> myDomPtLow, myDomPtUpper;
	PointI<TDim> myCenter;
	double myRadius {1.0};
	double minDistInitSegment {5.0};
	
	// Specific attributes to ImageMaskDomainCtr
	int myMaskThreshold {128};
	unsigned int myNbTry {100};
	double myMinDistanceToBorder {5.0};

public:
	Image myImage { Image(DomCT())} ;

	ImageD myDistanceImage { ImageD(DomCT()) };

	ImageMaskDomainCtrl(const ImageMaskDomainCtrl&)
	{
		std::cout << "copy domain!!" << std::endl;
	}

	ImageMaskDomainCtrl() {};
	
	// Constructor
	ImageMaskDomainCtrl(const std::string &fileImgDomain,
						int maskThreshold, PointI<TDim> ptRoot,
						unsigned int nbTry=100): myNbTry{nbTry}
	{
		myImage = DGtal::GenericReader<Image>::import(fileImgDomain,myMaskThreshold);
		myDistanceImage = GeomHelpers::getImageDistance<Image,ImageD>(myImage, myMaskThreshold);

		if ( !isInside(ptRoot) )
		{
			DGtal::trace.warning() << "ImageMaskDomainCtrl: Initial point given as input is not in domain." << std::endl;
			DGtal::trace.warning() << "ImageMaskDomainCtrl: Using default value from the maximal distant point." << std::endl;
			myCenter = maxDistantPointFromBorder();
		}
		else
		{
			myCenter = ptRoot;
		}
	}

	
	// Constructor
	ImageMaskDomainCtrl(const std::string &fileImgDomain,
						int maskThreshold, unsigned int nbTry=100)
	 : myNbTry(nbTry)
	{
		myImage = DGtal::GenericReader<Image>::import(fileImgDomain,myMaskThreshold);
		myDistanceImage = GeomHelpers::getImageDistance<Image,ImageD>(myImage, myMaskThreshold);
		myCenter = maxDistantPointFromBorder();
		checkImageDomain();
	};


	PointI<TDim>
	randomPoint()
	{
		bool found = false;
		PointI<TDim> pMin = myImage.domain().lowerBound();
		PointI<TDim> pMax = myImage.domain().upperBound();
		PointI<TDim> dp = pMax - pMin;
		PointI<TDim> pCand;
		unsigned int n = 0;

		while(!found && n < myNbTry)
		{
			for (unsigned int i = 0; i< TDim; i++)
			{
				pCand[i] = pMin[i]+(rand()%dp[i]);
			}

			found = isInside(pCand) && abs(myDistanceImage(pCand)) >= myMinDistanceToBorder;
			n++;
		}

		if (found)
		{
			return pCand;
		}

		if (n >= myNbTry)
		{
			for(auto p : myImage.domain())
			{
				if (isInside(p) && abs(myDistanceImage(p)) >= myMinDistanceToBorder)
				{
					return p;
				}
			}
		}
		return PointI<TDim>();
	}


	bool isInside(const PointI<TDim> &p) const
	{
		return myImage(p) > myMaskThreshold;
	}


	/**
	 * Check if the segment defined by two points intersect the domain.
	 *
	 * @param pt1 first point of the segment
	 * @param pt2  second point of the segment
	 */
	bool
	checkNoIntersectDomain(const PointI<TDim> &pt1, const PointI<TDim> &pt2) const
	{
		if (!myImage.domain().isInside(pt1)
			|| !myImage.domain().isInside(pt2))
		{
			return false;
		}

		PointD<TDim> dir = pt2 - pt1;
		dir /= dir.norm();
		PointD<TDim> p;

		for(unsigned int i=0; i<TDim; i++ ) { p[i]=pt1[i]; }

		for (unsigned int i = 0; i<(pt2 - pt1).norm(); i++)
		{
			PointD<TDim> p = pt1+dir*i;
			PointI<TDim> pI;
			for(unsigned int i=0; i<TDim; i++ ) { pI[i]=static_cast<int>(p[i]); }
			if (myImage(pI) < myMaskThreshold)
				return false;
		}

		return true;
	}
	

	PointI<TDim>
	maxDistantPointFromBorder() const 
	{
		double max = 0.0;
		PointI<TDim> pM;

		for(auto p: myDistanceImage.domain()) 
		{
			if (myDistanceImage(p) > max)
			{
				max = myDistanceImage(p);
				pM = p;
			}
		}

		return pM;
	}


	PointI<TDim>
	firstCandidatePoint() const 
	{
		PointI<TDim> res;
		bool find = searchRootFarthest(std::max(myDistanceImage(myCenter)/2.0, minDistInitSegment), res);
		assert(find);
		return res;
	}
	

	std::vector< std::vector< PointI<TDim> > >
	contours()
	{
		std::vector<std::vector< PointI<TDim> > > res;
		return res;
	}


	PointI<TDim>
	lowerBound() const
	{
		return myImage.domain().lowerBound();
	}


	PointI<TDim>
	upperBound() const
	{
		return myImage.domain().upperBound();
	}
	
private:
	// internal method
	
	bool
	searchRootFarthest(const double & d, PointI<TDim> &ptRoot ) const 
	{
		TDGset sPts = GeomHelpers::pointsOnSphere<PointI<TDim>, TDGset>(myCenter, d);

		for (const PointI<TDim> & p : sPts)
		{
			if (checkNoIntersectDomain(p, myCenter))
			{
				for(unsigned int i = 0; i < TDim; i++)
				{
					ptRoot[i] = p[i];
				}
				return true;
			}
		}
		return false;
	}
	

	void checkImageDomain()
	{
		bool isOk = false;
		// Check if at least one pixel of with foreground value exist:
		for (auto p: myImage.domain())
		{
			if (myImage(p) >= myMaskThreshold)
			{
				isOk = true;
				break;
			}
		}

		if (!isOk)
		{
			std::cout << "ImageMaskDomainCtrl: domain non valid since no point are inside the mask image." <<  std::endl;
		}
	}
};

template<>
std::vector<std::vector< PointI<2> > >
ImageMaskDomainCtrl<2>::contours()
{
	typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> TImage;
	typedef DGtal::functors::IntervalThresholder<typename TImage::Value> Binarizer;

	DGtal::Z2i::KSpace ks;

	if(! ks.init( myImage.domain().lowerBound(),
				 myImage.domain().upperBound(), true ))
	{
		DGtal::trace.error() << "Problem in KSpace initialisation"<< std::endl;
	}
	
	Binarizer b(myMaskThreshold, 255);
	DGtal::functors::PointFunctorPredicate<TImage,Binarizer> predicate(myImage, b);
	DGtal::trace.info() << "DGtal contour extraction from thresholds ["<<  myMaskThreshold << "," << 255 << "]" ;
	DGtal::SurfelAdjacency<2> sAdj( true );
	std::vector<std::vector< PointI<2> > > vectContoursBdryPointels;
	DGtal::Surfaces<DGtal::Z2i::KSpace>::extractAllPointContours4C( vectContoursBdryPointels, ks, predicate, sAdj );
	return vectContoursBdryPointels;
}


template<class DomCtr, int TDim>
std::vector< PointI<TDim> >
firstN_CandidatePoints(const DomCtr & domain_controller, unsigned int n)
{
	// create a distrubution of points at origin (points are placed around a n-D sphere)
	std::vector< PointD<TDim> > base_points(GeomHelpers::evenlySpreadPoints<TDim>(n));

	// use dichotomy to find the biggest sphere radius where all points fit inside the domain
	double r1 = domain_controller.myRadius;
	double r2 = 2 * r1;
	double precision = 0.5;

	// first we find r1 and r2 such that r2 = 2*r1 and all points with r1 are within the domain, and at least one point with r2 is not
	bool r_bounds_set = false;
	while(!r_bounds_set)
	{
		bool r1_points_in_dom = true;			// whether all points with r1 are inside the domain
		bool r2_points_in_dom = true;			// whether all points with r2 are inside the domain
		for(const PointD<TDim> p : base_points)
		{
			r1_points_in_dom = r1_points_in_dom 
				&& domain_controller.isInside(PointI<TDim>(domain_controller.myCenter + r1 * p));
			r2_points_in_dom = r2_points_in_dom 
				&& domain_controller.isInside(PointI<TDim>(domain_controller.myCenter + r2 * p));
		}

		r_bounds_set = r1_points_in_dom && !r2_points_in_dom;

		// adjust bounds if they don't satisfy criteria
		if(r1_points_in_dom)
		{	
			if(r2_points_in_dom)
			{
				r1 = r2;
				r2 *= 2;
			}
		}
		else
		{
			if(r2_points_in_dom)
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

	// dichotomy main loop
	while(fabs(r1 - r2) > precision)
	{
		double r_mean = (r1 + r2) / 2.0;

		double r_mean_within_dom = true;
		for(const PointD<TDim> & p : base_points)
		{
			r_mean_within_dom = r_mean_within_dom 
				&& domain_controller.isInside(PointI<TDim>(domain_controller.myCenter + r_mean * p));
		}

		// modify r1 or r2, while keeping their property
		if(r_mean_within_dom)
		{
			r1 = r_mean;
		}
		else
		{
			r2 = r_mean;
		}
	}

	std::vector< PointI<TDim> > res;

	for(const PointD<TDim> & p : base_points)
	{
		res.emplace_back(r1 * p);
	}

	return res;
}


///////////////////////////////////////////////////////////////////////////////

#endif // !defined DOMAIN_CONTROLLER_H

#undef DOMAIN_CONTROLLER_RECURSE
#endif // else defined(DOMAIN_CONTROLLER_RECURSE)