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
 *  -  TPoint  randomPoint():
 *    To get a point inside the domaine used in the tree construction.
 *
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

    CircularDomainCtrl(double radius, TPoint center)
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
};



/**
 *  Domain controller based on a Image Mask
 */
template<int TDim>
class ImageMaskDomainCtrl {
    typedef DGtal::PointVector<TDim, double> TPoint;
    typedef  DGtal::PointVector<TDim, int> TPointI;

    typedef DGtal::SpaceND< TDim, int >   SpaceCT;
    typedef DGtal::HyperRectDomain<SpaceCT> DomCT;
    typedef DGtal::ImageContainerBySTLVector< DomCT, int> Image;

    TPoint myDomPtLow, myDomPtUpper;
    int myMaskThreshold {128};
    unsigned int myNbTry {0};
    Image myImage;
    Image myDistanceImage;
    // Constructor for Masked domain type
    ImageMaskDomainCtrl(const std::string &fileImgDomain,
                        int maskThreshold, unsigned int nbTry=100): myNbTry{nbTry}{
        myImage = DGtal::GenericReader<Image>::import( fileImgDomain );
        myDistanceImage = GeomHelpers::getImageDistance<Image, Image>(myImage);
    };
    
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
    
    
};



///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods 
#include "DomainController.ih"



///////////////////////////////////////////////////////////////////////////////

#endif // !defined DOMAIN_CONTROLLER_H

#undef DOMAIN_CONTROLLER_RECURSE
#endif // else defined(DOMAIN_CONTROLLER_RECURSE)





