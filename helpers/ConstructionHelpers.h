
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
#include "DGtal/io/readers/GenericReader.h"
#include "CoronaryArteryTree.h"



namespace ConstructionHelpers {



/**
 * Helpers fonction to construction the tree. (mainly from circular domain and image organ for testing).
 */
template<int TDim>
inline
void constructTree(double aPerf, int nbTerm,
                   std::string imageOrgan, unsigned int fgTh = 128,
                   bool verbose = false, DGtal::PointVector<TDim, int> ptCenter = DGtal::PointVector<TDim, int>::diagonal(0),
                   unsigned int distSearchRoot = 10){
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test random adds with distance constraint");
  srand (time(NULL));
  double rRoot = 1.0;
  std::string filename;
  
  CoronaryArteryTree<TDim> cTree (aPerf, nbTerm, rRoot, ptCenter);
  if (imageOrgan != ""){
    auto img = DGtal::GenericReader<typename CoronaryArteryTree<TDim>::Image>::import( imageOrgan );
    bool restrainedOK = cTree.restrainDomain(img, fgTh);
    if (restrainedOK){
      DGtal::trace.info() << "Using restrained image  " << imageOrgan << std::endl;
      cTree.myVectSegments[1].myCoordinate = ptCenter;
      cTree.myTreeCenter = ptCenter;
      if (ptCenter[0]==-1 && ptCenter[0]==-1) {
        auto pC = cTree.getDomainCenter();
        cTree.myVectSegments[1].myCoordinate = pC;
      }
      DGtal::PointVector<TDim, double> pRoot;
      cTree.searchRootFarthest(distSearchRoot, pRoot);
      cTree.myVectSegments[0].myCoordinate = pRoot;
    }
  }
  bool isOK = false;
  unsigned int nbSeed = cTree.my_NTerm;
  for (unsigned int i = 1; i < nbSeed; i++) {
    DGtal::trace.progressBar(i, nbSeed);
    int nbSol = 0, itOpt = 0;
    CoronaryArteryTree<TDim> cTreeOpt = cTree;
    double volOpt = -1.0, vol = 0.0;
    while (nbSol==0) {
      auto pt = cTree.generateNewLocation(100);
      std::vector<unsigned int> vecN = cTree.getN_NearestSegments(pt,cTree.myNumNeighbor);
      for(size_t it=0; it<vecN.size(); it++) {
        //if(!cTree.isIntersecting(pt, cTree.findBarycenter(pt, vecN.at(it)),vecN.at(it),n))
        if(!cTree.isIntersecting(pt, cTree.findBarycenter(pt, vecN.at(it)),vecN.at(it),cTree.myNumNeighbor, 2*cTree.myVectSegments[vecN.at(it)].myRadius)) {
          CoronaryArteryTree<TDim> cTree1 = cTree;
          isOK = cTree1.isAddable(pt,vecN.at(it), 100, 0.01, cTree1.myNumNeighbor, verbose);
          if(isOK) {
            vol = cTree1.computeTotalVolume(1);
            if(volOpt<0.0) {
              volOpt = vol;
              cTreeOpt = cTree1;
              itOpt = it;
            }
            else {
              if(volOpt>vol) {
                volOpt = vol;
                cTreeOpt = cTree1;
                itOpt = it;
              }
            }
            nbSol++;
          }
        }
      }
    }
    cTree = cTreeOpt;
    cTree.updateLengthFactor();
    cTree.updateResistanceFromRoot();
    cTree.updateRootRadius();
  }
  if (verbose){
    std::cout<<"====> Aperf="<<cTree.myRsupp*cTree.myRsupp*cTree.my_NTerm*M_PI<<" == "<<aPerf<<std::endl;
  }
  filename = "testCCO_"+std::to_string(nbTerm)+".eps";
  cTree.exportBoardDisplay(filename.c_str(), 1.0);
  cTree.exportBoardDisplay("result.eps", 1.0);
  cTree.exportBoardDisplay("result.svg", 1.0);

  cTree.myBoard.clear();
}

/**
 * Helpers fonction to construction the tree with autoSearch of the center and root.
 * The center is defined from the maximal distance map and the root point is searched on the image domain.
 */
template<typename TPoint>
inline
void constructTreeImageDomain2D(double aPerf, int nbTerm,
                              std::string imageOrgan, unsigned int fgTh = 128,
                              bool verbose = false){
  // searching center from maximak distance.
  auto img = DGtal::GenericReader<typename CoronaryArteryTree<TPoint::dimension>::Image>::import( imageOrgan );
  auto imgDist = GeomHelpers::getImageDistance2D<typename CoronaryArteryTree<TPoint::dimension>::Image,
                                               typename CoronaryArteryTree<TPoint::dimension>::ImageDist>(img);
  double m = 0.0;
  DGtal::PointVector<2, int> pM;
  for(auto p: imgDist.domain()) {if (imgDist(p) > m ){m = imgDist(p); pM =  DGtal::PointVector<2, int>(p[0], p[1]);}}
  if (verbose){
    DGtal::trace.info() << "center point found: " << pM << "maximal value:"
                        << m <<   std::endl;
  }
  constructTree<2>(aPerf, nbTerm, imageOrgan, fgTh, verbose, pM,
                static_cast<unsigned int >(m/2.0));
}


/**
 * Helpers fonction to construction the tree with autoSearch of the center and root.
 * The center is defined from the maximal distance map and the root point is searched on the image domain.
 */
template<typename TPoint>
inline
void constructTreeImageDomain3D(double aPerf, int nbTerm,
                              std::string imageOrgan, unsigned int fgTh = 128,
                              bool verbose = false){
  // searching center from maximak distance.
  auto img = DGtal::GenericReader<typename CoronaryArteryTree<TPoint::dimension>::Image>::import( imageOrgan );
  auto imgDist = GeomHelpers::getImageDistance3D<typename CoronaryArteryTree<TPoint::dimension>::Image,
                                               typename CoronaryArteryTree<TPoint::dimension>::ImageDist>(img);
  double m = 0.0;
  TPoint pM;
  for(auto p: imgDist.domain()) {if (imgDist(p) > m ){m = imgDist(p); pM = TPoint(p[0], p[1]);}}
  if (verbose){
    DGtal::trace.info() << "center point found: " << pM << "maximal value:"
                        << m <<   std::endl;
  }
  constructTree<3>(aPerf, nbTerm, imageOrgan, fgTh, verbose, pM,
                static_cast<unsigned int >(m)/2.0);
}

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




