
#pragma once

#if defined(CONSTRUCTIONHELPERS_RECURSES)
#error Recursive header files inclusion detected in ConstructionHelpers.h
#else // defined(CONSTRUCTIONHELPERS_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CONSTRUCTIONHELPERS_RECURSES

#if !defined CONSTRUCTIONHELPERS_h
/** Prevents repeated inclusion of headers. */
#define CONSTRUCTIONHELPERS_h

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/kernel/BasicPointPredicates.h"


namespace ConstructionHelpers {

/**
 * Helpers fonction to construction the tree with autoSearch of the center and root.
 * The center is defined from the maximal distance map and the root point is searched on the image domain.
 */
void constructTreeImageDomain(double aPerf, int nbTerm,
                   std::string imageOrgan, unsigned int fgTh = 128,
                   bool verbose = false);



/**
 * Helpers fonction to construction the tree. (mainly from circular domain and image organ for testing).
 */
void constructTree(double aPerf, int nbTerm,
                   std::string imageOrgan, unsigned int fgTh = 128,
                   bool verbose = false, DGtal::Z2i::RealPoint ptCenter= DGtal::Z2i::RealPoint(0,0),
                   unsigned int distSearchRoot = 10);


template< typename TImage>
inline
std::vector<std::vector<DGtal::Z2i::Point> > getImageContours(const TImage &image,
                                                              unsigned int threshold=128){
  typedef DGtal::functors::IntervalThresholder<typename TImage::Value> Binarizer;
  DGtal::Z2i::KSpace ks;
  if(! ks.init( image.domain().lowerBound(),
               image.domain().upperBound(), true )){
    DGtal::trace.error() << "Problem in KSpace initialisation"<< std::endl;
  }
  
  Binarizer b(threshold, 255);
  DGtal::functors::PointFunctorPredicate<TImage,Binarizer> predicate(image, b);
  DGtal::trace.info() << "DGtal contour extraction from thresholds ["<<  threshold << "," << 255 << "]" ;
  DGtal::SurfelAdjacency<2> sAdj( true );
  std::vector< std::vector< DGtal::Z2i::Point >  >  vectContoursBdryPointels;
  DGtal::Surfaces<DGtal::Z2i::KSpace>::extractAllPointContours4C( vectContoursBdryPointels, ks, predicate, sAdj );
  return vectContoursBdryPointels;
}


}

#endif // !defined CONSTRUCTIONHELPERS_h

#undef CONSTRUCTIONHELPERS_RECURSES
#endif // else defined(CONSTRUCTIONHELPERS_RECURSES)




