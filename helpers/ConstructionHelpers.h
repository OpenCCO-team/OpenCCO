
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

void constructTree(double aPerf, int nbTerm,
                   std::string imageOrgan, unsigned int fgTh = 128, int minDBorder = 0,
                   bool verbose = false, DGtal::Z2i::Point ptCenter= DGtal::Z2i::Point(0,0));

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




