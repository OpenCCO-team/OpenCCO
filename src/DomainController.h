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
 *  - myCenter: used to orient the construction according the domain.
 *  - bool isInside(const TPoint &p):
 *    Test if a point is inside the domain or not (a point can be outside after the optmisation).
 *  - TPoint  randomPoint():
 *    To get a point inside the domaine used in the tree construction.
 *  - bool checkNoIntersectDomain(const TPoint &pt1, const TPoint &pt2):
 Usefull to ensure is a whole segment is inside the domain.
 *  - maxDistantPointFromBorder() const:
 *    usefull to determine a starting point and help to construct the tree.
 *  - TPoint firtCandidatePoint() const :
 *    use to initiate the reconstruction of the tree.
 **/


/**
 *  Domain controller based on a ball
 */
template<int TDim>
class CircularDomainCtrl {
public:
    typedef DGtal::PointVector<TDim, double> TPoint;
    double myRadius {1.0};
    TPoint myCenter;
    
    // Fixme to be removed since no need when finalized
    bool myIsImageDomainRestrained = false;
    // Constructor for ImplicitCirc type
    CircularDomainCtrl(){};
    
    CircularDomainCtrl(double radius, const TPoint &center)
    {
        myRadius = radius;
        myCenter = center;
    };
    
    bool isInside(const TPoint &p)
    {
        return (myCenter-p).norm() < myRadius;
    }
    
    TPoint
    randomPoint() {
        bool found = false;
        TPoint p;
        while(!found){
            double ss = 0.0;
            for(unsigned int i = 0;i<TDim; i++ ){
                p[i] = ((double)rand() / RAND_MAX)*2.0*myRadius - myRadius;
            }
            found = isInside(p);
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
    
    TPoint
    firtCandidatePoint() const {
        TPoint res;
        if (TDim == 2){
            res[0] = 0.0;
            res[1] = myRadius;
        }else  if (TDim == 3){
            res[0] = 0.0;
            res[1] = myRadius;
            res[2] = 0.0;
            
        }
        
        return res;
        
    }
    std::vector<std::vector< TPoint > >
    contours()
    {
        std::vector<std::vector< TPoint > > res;
        return res;
    }
    TPoint
    lowerBound()
    {
        TPoint p = TPoint::diagonal(myRadius*0.01);
        return myCenter - p;
    }
    TPoint
    upperBound()
    {
        TPoint p = TPoint::diagonal(myRadius*0.01);
        return myCenter + p;
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
    typedef typename DGtal::DigitalSetSelector<DomCT,
    DGtal::BIG_DS+
    DGtal::HIGH_BEL_DS>::Type TDGset;
    
    
    TPoint myDomPtLow, myDomPtUpper;
    int myMaskThreshold {128};
    unsigned int myNbTry {100};
    // Fixme to be removed since no need when finalized
    bool myIsImageDomainRestrained = true;
    TPointI myCenter;
    double minDistInitSegment {5.0};

    
public:
    Image myImage {Image(DomCT())} ;
    ImageD myDistanceImage {ImageD(DomCT())};
    
    ImageMaskDomainCtrl(){};
    
    // Constructor
    ImageMaskDomainCtrl(const std::string &fileImgDomain,
                        int maskThreshold, TPointI ptRoot,
                        unsigned int nbTry=100): myNbTry{nbTry}
    {
        myImage = DGtal::GenericReader<Image>::import(fileImgDomain,myMaskThreshold);
        myDistanceImage = GeomHelpers::getImageDistance<Image,ImageD>(myImage,
                                                                      myMaskThreshold );
        if ( !isInside(ptRoot) ){
            DGtal::trace.warning() << "ImageMaskDomainCtrl: Initial point given as input is not in domain." << std::endl;
            DGtal::trace.warning() << "ImageMaskDomainCtrl: Using default value from the maximal distant point." << std::endl;
            myCenter = maxDistantPointFromBorder();
        }else{
            myCenter = ptRoot;
        }

    }
    
    // Constructor
    ImageMaskDomainCtrl(const std::string &fileImgDomain,
                        int maskThreshold, unsigned int nbTry=100):
                                                    myNbTry{nbTry}
    {
        myImage = DGtal::GenericReader<Image>::import(fileImgDomain,myMaskThreshold);
        myDistanceImage = GeomHelpers::getImageDistance<Image,ImageD>(myImage,
                                                                      myMaskThreshold );
        myCenter = maxDistantPointFromBorder();
        checkImageDomain();
    };
    
    
    
    
    
    TPointI
    randomPoint()
    {
        bool found = false;
        TPointI pMin = myImage.domain().lowerBound();
        TPointI pMax = myImage.domain().upperBound();
        TPointI dp = pMax - pMin;
        TPointI pCand;
        unsigned int n = 0;
        while(!found && n < myNbTry){
            for (unsigned int i = 0; i< TDim; i++){
                pCand[i] = pMin[i]+(rand()%dp[i]);
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
    bool isInside(const TPointI &p){
        return myImage(p) > myMaskThreshold;
    }
    /**
     * Check if the segment defined by two points intersect the domain.
     *
     * @param pt1 first point of the segment
     * @param pt2  second point of the segment
     */
    bool
    checkNoIntersectDomain(const TPointI &pt1, const TPointI &pt2) const
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
            TPointI pI;
            for(unsigned int i=0; i<TPoint::dimension; i++ ){pI[i]=static_cast<int>(p[i]);}
            if (myImage(pI) < myMaskThreshold)
                return false;
        }
        return true;
    }
    
    TPointI
    maxDistantPointFromBorder() const {
        double m = 0.0;
        TPointI pM;
        for(auto p: myDistanceImage.domain()) {
            if (myDistanceImage(p) > m ){m = myDistanceImage(p);
                for (unsigned int i=0; i<TDim; i++){
                    pM[i] = p[i];
                }
            }
        }
        return pM;
    }
    TPoint
    firtCandidatePoint() const {
        TPointI res;
        searchRootFarthest(std::max(myDistanceImage(myCenter)/2.0, minDistInitSegment), res);
        return res;
    }
    
    std::vector<std::vector< TPointI > >
    contours()
    {
        std::vector<std::vector< TPointI > > res;
        return res;
    }
    TPointI
    lowerBound()
    {
        
        return myImage.domain().lowerBound();
    }
    TPointI
    upperBound()
    {
        return myImage.domain().upperBound();
    }
    
private:
    // internal method
    
    bool
    searchRootFarthest(const double & d, TPointI &ptRoot ) const {
        typedef DGtal::SpaceND<TDim, int> Space;
        //        auto sPts =  GeomHelpers::pointsOnSphere(ptRoot, d);
        DomCT aDom;
        TDGset sPts = GeomHelpers::pointsOnSphere<TPointI, TDGset>(myCenter, d);
        for (const TPointI &p : sPts){
            if (checkNoIntersectDomain(p, myCenter)){
                for(unsigned int i = 0; i < TDim; i++){
                    ptRoot[i] = p[i];
                }
                return true;
            }
        }
        return false;
    }
    
    void checkImageDomain(){
        bool isOk = false;
        // Check if at least one pixel of with foreground value exist:
        for (auto p: myImage.domain()){
            if (myImage(p) >= myMaskThreshold){
                isOk = true;
                break;
            }
        }
        if (!isOk){
            std::cout << "ImageMaskDomainCtrl: domain non valid since no point are inside the mask image." <<  std::endl;
        }
    }
    
};

template<>
std::vector<std::vector< typename ImageMaskDomainCtrl<2>::TPointI > >
ImageMaskDomainCtrl<2>::contours()
{
    typedef  DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> TImage;
    typedef DGtal::functors::IntervalThresholder<typename TImage::Value> Binarizer;
    DGtal::Z2i::KSpace ks;
    if(! ks.init( myImage.domain().lowerBound(),
                 myImage.domain().upperBound(), true )){
        DGtal::trace.error() << "Problem in KSpace initialisation"<< std::endl;
    }
    
    Binarizer b(myMaskThreshold, 255);
    DGtal::functors::PointFunctorPredicate<TImage,Binarizer> predicate(myImage, b);
    DGtal::trace.info() << "DGtal contour extraction from thresholds ["<<  myMaskThreshold << "," << 255 << "]" ;
    DGtal::SurfelAdjacency<2> sAdj( true );
    std::vector<std::vector< typename ImageMaskDomainCtrl<2>::TPointI > > vectContoursBdryPointels;
    DGtal::Surfaces<DGtal::Z2i::KSpace>::extractAllPointContours4C( vectContoursBdryPointels, ks, predicate, sAdj );
    return vectContoursBdryPointels;
}



///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods 
#include "DomainController.ih"



///////////////////////////////////////////////////////////////////////////////

#endif // !defined DOMAIN_CONTROLLER_H

#undef DOMAIN_CONTROLLER_RECURSE
#endif // else defined(DOMAIN_CONTROLLER_RECURSE)





