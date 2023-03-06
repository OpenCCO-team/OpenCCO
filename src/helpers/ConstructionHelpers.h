
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


template<typename DomCtr, int TDim>
inline
void
initFirtElemTree(CoronaryArteryTree< DomCtr, TDim > &aTree,
                 unsigned int distSearchRoot, bool verbose = false){
    DGtal::PointVector<TDim, double> pRoot;
    bool restrainedOK = aTree.restrainDomain(aTree.myDomainController.myImage,              aTree.myForegroundThreshold);
    if (restrainedOK){
        DGtal::trace.info() << "Using restrained image  "  << std::endl;
    }
    aTree.searchRootFarthest(distSearchRoot, pRoot);
    aTree.myVectSegments[0].myCoordinate = pRoot;
}


template<typename DomCtr, int TDim>
inline
void
expandTree(CoronaryArteryTree< DomCtr, TDim > &aTree, bool verbose = false)
{
    srand ((unsigned int) time(NULL));
    
    bool isOK = false;
     unsigned int nbSeed = aTree.my_NTerm;
     for (unsigned int i = 1; i < nbSeed; i++) {
       DGtal::trace.progressBar(i, nbSeed);
       size_t nbSol = 0, itOpt = 0;
       CoronaryArteryTree< DomCtr, TDim > cTreeOpt = aTree;
       double volOpt = -1.0, vol = 0.0;
       while (nbSol==0) {
         auto pt = aTree.generateNewLocation(100);
         std::vector<unsigned int> vecN = aTree.getN_NearestSegments(pt,aTree.myNumNeighbor);
         for(size_t it=0; it<vecN.size(); it++) {
           if(!aTree.isIntersecting(pt, aTree.findBarycenter(pt, vecN.at(it)),vecN.at(it),aTree.myNumNeighbor, 2*aTree.myVectSegments[vecN.at(it)].myRadius)) {
             CoronaryArteryTree< DomCtr, TDim  > cTree1 = aTree;
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
         aTree = cTreeOpt;
         aTree.updateLengthFactor();
         aTree.updateResistanceFromRoot();
         aTree.updateRootRadius();
     }
     if (verbose){
       std::cout<<"====> Aperf="<<aTree.myRsupp*aTree.myRsupp*aTree.my_NTerm*M_PI<<" == "<<aTree.my_aPerf<<std::endl;
     }
}



/**
 * Helpers fonction to construction the tree with autoSearch of the center and root.
 * The center is defined from the maximal distance map and the root point is searched on the image domain.
 */
template<typename TPoint, typename DomCtr>
inline
CoronaryArteryTree<DomCtr, 2> constructTreeImageDomain2D(double aPerf, int nbTerm,
                              std::string imageOrgan, unsigned int fgTh = 128,
                              bool verbose = false){
  // searching center from maximak distance.
  auto img = DGtal::GenericReader<typename CoronaryArteryTree<DomCtr,TPoint::dimension>::Image>::import( imageOrgan );
  auto imgDist = GeomHelpers::getImageDistance<typename CoronaryArteryTree<DomCtr, TPoint::dimension>::Image,
                                               typename CoronaryArteryTree<DomCtr, TPoint::dimension>::ImageDist>(img);
  double m = 0.0;
  DGtal::PointVector<2, int> pM;
  for(auto p: imgDist.domain()) {if (imgDist(p) > m ){m = imgDist(p); pM =  DGtal::PointVector<2, int>(p[0], p[1]);}}
  if (verbose){
    DGtal::trace.info() << "center point found: " << pM << "maximal value:"
                        << m <<   std::endl;
  }
    return constructTree<DomCtr, 2>(aPerf, nbTerm, imageOrgan, fgTh, verbose, pM,
                static_cast<unsigned int >(m/2.0));
}


/**
 * Helpers fonction to construction the tree with autoSearch of the center and root.
 * The center is defined from the maximal distance map and the root point is searched on the image domain.
 */
template<typename TPoint, typename DomCtr>
inline
CoronaryArteryTree<DomCtr, 3> constructTreeImageDomain3D(double aPerf, int nbTerm,
                              std::string imageOrgan, unsigned int fgTh = 128,
                              bool verbose = false){
  // searching center from maximak distance.
  auto img = DGtal::GenericReader<typename CoronaryArteryTree<DomCtr, TPoint::dimension>::Image>::import( imageOrgan );
  auto imgDist = GeomHelpers::getImageDistance<typename CoronaryArteryTree<DomCtr,3>::Image,
                                               typename CoronaryArteryTree<DomCtr, 3>::ImageDist>(img,fgTh);
  double m = 0.0;
  DGtal::PointVector<3, int> pM;
  for(auto p: imgDist.domain()) {
    if (imgDist(p) > m ){m = imgDist(p); pM = DGtal::PointVector<3, int>(p[0], p[1], p[2]);}}
  if (verbose){
    DGtal::trace.info() << "center point found: " << pM << "maximal value:"
                        << m <<   std::endl;
  }
  return constructTree<DomCtr, 3>(aPerf, nbTerm, imageOrgan, fgTh, verbose, pM,
                static_cast<unsigned int >(m/2.0));
}


template<typename TImage >
inline
std::vector<std::vector<typename TImage::Domain::Point > >
getImageContours(const TImage &image,
                 unsigned int threshold=128){
  std::vector<std::vector<typename TImage::Domain::Point > > v;
  DGtal::trace.error() << "Use CCO is only implemented in 2D and 3D and use ImageContainerBySTLVector with unsigned char"
                << "to export domain."
                << "You use such an image : " << image <<  std::endl;
      throw 1;
  return v;
}

/**
  Template Specialisation in 2D  to export image contour of the restricted image domain.
 */
template< >
inline
std::vector<std::vector<DGtal::Z2i::Point > >
getImageContours(const DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> &image,
                 unsigned int threshold){
    typedef  DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> TImage;
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




/**
  Template Specialisation in 3D  to export surfel image  border of the restricted image domain.
 * Todo @BK
 */
template< >
inline
std::vector<std::vector<DGtal::Z3i::Point > >
getImageContours(const DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, unsigned char> &image,
                 unsigned int threshold){
  std::vector<std::vector<DGtal::Z3i::Point > > v;
  
  return v;
}




}

#endif // !defined CONSTRUCTIONHELPERS_h

#undef CONSTRUCTIONHELPERS_RECURSES
#endif // else defined(CONSTRUCTIONHELPERS_RECURSES)




