
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
    typedef typename CoronaryArteryTree< DomCtr, TDim >::TreeState TState;
    srand ((unsigned int) time(NULL));
    bool isOK = false;
     unsigned int nbSeed = aTree.bParam.my_NTerm;
    for (unsigned int i = 1; i < nbSeed; i++) {
        DGtal::trace.progressBar(i, nbSeed);
        size_t nbSol = 0, itOpt = 0;
        TState stateOpt = aTree.state();
        TState stateBase = aTree.state();

        double volOpt = -1.0, vol = 0.0;
        while (nbSol==0) {
            auto pt = aTree.generateNewLocation(100);
            std::vector<unsigned int> vecN = aTree.getN_NearestSegments(pt,aTree.iParam.myNumNeighbor);
            for(size_t it=0; it<vecN.size(); it++) {
                if(!aTree.isIntersecting(pt, aTree.findBarycenter(pt, vecN.at(it)),vecN.at(it),aTree.iParam.myNumNeighbor, 2*aTree.myVectSegments[vecN.at(it)].myRadius)) {
                    // CTree1 copie d'expé
                    aTree.restaureState(stateBase);
                    isOK = aTree.isAddable(pt,vecN.at(it), 100, 0.01, aTree.iParam.myNumNeighbor, verbose);
                    if(isOK) {
                        vol = aTree.computeTotalVolume(1);
                        if(volOpt<0.0) {
                            volOpt = vol;
                            stateOpt = aTree.state();
                            itOpt = it;
                        }
                        else {
                            if(volOpt>vol) {
                                volOpt = vol;
                                stateOpt = aTree.state();
                                itOpt = it;
                            }
                        }
                        nbSol++;
                    }
                }
            }
        }
        aTree.restaureState(stateOpt);
        aTree.updateLengthFactor();
        aTree.updateResistanceFromRoot();
        aTree.updateRootRadius();
    }
     if (verbose){
       std::cout<<"====> Aperf="<<aTree.iParam.myRsupp*aTree.iParam.myRsupp*aTree.bParam.my_NTerm*M_PI<<" == "<<aTree.bParam.my_aPerf<<std::endl;
     }
}






template<typename TImage >
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




