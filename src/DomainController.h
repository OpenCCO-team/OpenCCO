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

/**
 * Classes to control the construction domain of the coronary aretery tree.
 *  It is used during the reconstruction using the following points:
 *  - bool isInside(const TPoint &p):
 *    Test if a point is inside the domain or not (a point can be outside after the optmisation).
 *  - TPoint  randomPoint():
 *    To get a point inside the domaine used in the tree construction.
 *  - bool checkNoIntersectDomain(const TPoint &pt1, const TPoint &pt2):
      Usefull to ensure is a whole segment is inside the domain.
 *  - maxDistantPointFromBorder() const:
 *    usefull to determine a starting point and help to construct the tree.

 */


/**
 *  Domain controller based on a ball
 */
template<int TDim>
class CircularDomainCtrl {
public:
    typedef DGtal::PointVector<TDim, double> TPoint;
    double myRadius {0.0};
    TPoint myCenter;
    // Constructor for ImplicitCirc type
    CircularDomainCtrl(){};
    
    CircularDomainCtrl(double radius, const TPoint &center)
    {
        myRadius = radius;
        myCenter = center;
    };
    
    bool isInside(const TPoint &p)
    {
        return p.norm() < myRadius;
    }
    
    TPoint
    randomPoint() {
      bool found = false;
      TPoint p;
      while(!found){
          double ss = 0.0;
          for(unsigned int i = 0;i<TDim; i++ ){
              p[i] = ((double)rand() / RAND_MAX)*2.0*myRadius - myRadius;
              ss += p[i]*p[i];
          }
        found = ss < myRadius*myRadius;
      }
      return p + myCenter;
    }
    /**
     * Get the supported domain of the tree. By default it is defined from the circle center.
     * If the domain if defined from a mask image, the center if computed from the imate center.
     */
    TPoint getDomainCenter() const{
        return myCenter;
    }
    bool
    checkNoIntersectDomain(const TPoint &pt1, const TPoint &pt2)
    {
        return isInside(pt1) && isInside(pt2);
    }
    
    TPoint
    maxDistantPointFromBorder() const {
        return myCenter;
    }
};



/**
 *  Domain controller based on a Image Mask
 */
template<int TDim>
class ImageMaskDomainCtrl {
public:
    typedef DGtal::PointVector<TDim, double> TPoint;
    typedef  DGtal::PointVector<TDim, int> TPointI;

    typedef DGtal::SpaceND< TDim, int >   SpaceCT;
    typedef DGtal::HyperRectDomain<SpaceCT> DomCT;
    typedef DGtal::ImageContainerBySTLVector< DomCT, unsigned char> Image;
    typedef DGtal::ImageContainerBySTLVector< DomCT, double> ImageD;

    
    TPoint myDomPtLow, myDomPtUpper;
    int myMaskThreshold {128};
    unsigned int myNbTry {0};
   
public:
    Image myImage;
    ImageD myDistanceImage;
    
    ImageMaskDomainCtrl(): myImage{Image(DomCT())},
                           myDistanceImage {ImageD(DomCT())} {};
    
    
    // Constructor for Masked domain type
    ImageMaskDomainCtrl(const std::string &fileImgDomain,
                        int maskThreshold, unsigned int nbTry=100): myNbTry{nbTry},       myImage {DGtal::GenericReader<Image>::import( fileImgDomain,myMaskThreshold )},
    myDistanceImage { GeomHelpers::getImageDistance3D<Image,ImageD>(myImage,myMaskThreshold )}
    {};
    
    TPointI
    randomPoint()
    {
      bool found = false;
      unsigned int x = 0;
      unsigned int y = 0;
      TPointI pMin = myImage.domain().lowerBound();
      TPointI pMax = myImage.domain().upperBound();
      TPointI dp = pMax - pMin;
      TPointI pCand;
      unsigned int n = 0;
      while(!found && n < myNbTry){
          for (unsigned int i = 0; i< TDim; i++){
              pCand[i] = pCand[i]%dp[i];
          }
        found = myImage(pCand)>=myMaskThreshold &&
                abs(myDistanceImage(pCand)) >= 10.0;
        n++;
      }
      if (n >= myNbTry){
        for(auto p : myImage.domain()){if (myImage(p)>=myMaskThreshold && abs(myDistanceImage(p)) >= 10.0 ) return p;}
      }
      return pCand;
    }
    bool isInside(const TPoint &p){
        return myImage(p) > myMaskThreshold;
    }
    /**
     * Check if the segment defined by two points intersect the domain.
     *
     * @param pt1 first point of the segment
     * @param pt2  second point of the segment
     */
   bool
   checkNoIntersectDomain(const TPoint &pt1, const TPoint &pt2)
    {
      if (!myImage.domain().isInside(pt1) ||
          !myImage.domain().isInside(pt2)){
        return false;
      }
      DGtal::PointVector<TPoint::dimension, double> dir = pt2 - pt1;
      dir /= dir.norm();
      DGtal::PointVector<TPoint::dimension, double>  p;
      for(unsigned int i=0; i<TPoint::dimension; i++ ){p[i]=pt1[i];}
      for (unsigned int i = 0; i<(pt2 - pt1).norm(); i++){
        DGtal::PointVector<TPoint::dimension, double>  p = pt1+dir*i;
        TPoint pI;
        for(unsigned int i=0; i<TPoint::dimension; i++ ){pI[i]=static_cast<int>(p[i]);}
        if (myImage(pI) < myMaskThreshold)
          return false;
      }
      return true;
    }
    
    TPointI
    maxDistantPointFromBorder() const {
        double m = 0.0;
        DGtal::PointVector<3, int> pM;
        for(auto p: myDistanceImage.domain()) {
            if (myDistanceImage(p) > m ){m = myDistanceImage(p);
                pM = DGtal::PointVector<3, int>(p[0], p[1], p[2]);}
            
        }
        return pM;
    }
    
};



///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods 
#include "DomainController.ih"



///////////////////////////////////////////////////////////////////////////////

#endif // !defined DOMAIN_CONTROLLER_H

#undef DOMAIN_CONTROLLER_RECURSE
#endif // else defined(DOMAIN_CONTROLLER_RECURSE)





