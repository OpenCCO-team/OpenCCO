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




template<int TDim>
class DomainCtrl {
    // point is let in double to simplify vectorial domain and template process.
    typedef DGtal::PointVector<TDim, double> TPoint;
    
    virtual bool isInside(const TPoint &p);
    virtual TPoint randomPoint();
};


/**
    Domain controller based on a ball
 */

template<int TDim>
class CircularDomainCtrl: DomainCtrl<TDim> {
public:
    typedef  typename DomainCtrl<TDim>::TPoint TPoint;
    double myRadius {0.0};
    TPoint myCenter;
    // Constructor for ImplicitCirc type

    CircularDomainCtrl(double radius, TPoint center){
        myRadius = radius;
        myCenter = center;
    };
    
    bool isInside(const TPoint &p){
        return p.norm() < myRadius;
    }
    
    TPoint
    randomPoint()
    {
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
        return p + myCenter;;
    }
};



/**
    Domain controller based on a Image Mask
 */
template<int TDim>
class ImageMaskDomainCtrl: DomainCtrl<TDim> {
    typedef  typename DomainCtrl<TDim>::TPoint TPoint;
    typedef  DGtal::PointVector<TDim, int> TPointI;

    typedef DGtal::SpaceND< TDim, int >   SpaceCT;
    typedef DGtal::HyperRectDomain<SpaceCT> DomCT;
    typedef DGtal::ImageContainerBySTLVector< DomCT, int> Image;

    TPoint myDomPtLow, myDomPtUpper;
    int myMaskThreshold {128};
    Image myImage;
    Image myDistanceImage;
    // Constructor for Masked domain type
    ImageMaskDomainCtrl(const std::string &fileImgDomain, int maskThreshold){
        myImage = DGtal::GenericReader<Image>::import( fileImgDomain );
        myDistanceImage = GeomHelpers::getImageDistance<Image, Image>(myImage);
    };
    
    bool isInside(const TPoint &p){
        
        return myImage(p) > myMaskThreshold;
    }
    
    TPoint
    randomPoint(unsigned int nbTry=0){
        bool found = false;
        unsigned int x = 0;
        unsigned int y = 0;
        unsigned int z = 0;
        auto pMin = myImage.domain().lowerBound();
        auto pMax = myImage.domain().upperBound();
        auto dp = pMax - pMin;

        DGtal::PointVector<TDim, int > pCand;
        unsigned int n = 0;
        while(!found && n < nbTry){
            TPoint p;
            for (unsigned int i = 0; i < TDim; i++ )
            {
                p[i] = rand()%dp[i];
                p[i+1] = rand()%dp[i+1];
                p[i+2] = rand()%dp[i+2];
            }
          pCand[0] = static_cast<int>(pMin[0] +x);
          pCand[1] = static_cast<int>(pMin[1] +y);
          pCand[2] = static_cast<int>(pMin[2] +z);
          found = myImage(pCand)>= myMaskThreshold  && abs(myDistanceImage(pCand)) >= 10.0;
          n++;
        }
        if (n >= nbTry){
          for(auto p : myImage.domain()){
              if (myImage(p)>=myMaskThreshold &&
                  abs(myDistanceImage(p)) >= 10.0 )
                  return p;
          }
        }
        return pCand;
    }
};



///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods 
#include "DomainController.ih"



///////////////////////////////////////////////////////////////////////////////

#endif // !defined DOMAIN_CONTROLLER_H

#undef DOMAIN_CONTROLLER_RECURSE
#endif // else defined(DOMAIN_CONTROLLER_RECURSE)





