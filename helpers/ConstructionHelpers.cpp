
#include "ConstructionHelpers.h"
#include "CoronaryArteryTree.h"
#include "DGtal/geometry/helpers/ContourHelper.h"
#include "DGtal/io/readers/GenericReader.h"

void
ConstructionHelpers::constructTreeImageDomain(double aPerf, int nbTerm,
                                              std::string imageOrgan,
                                              unsigned int fgTh,
                                              bool verbose){
  // searching center from maximak distance.
  auto img = DGtal::GenericReader<CoronaryArteryTree::Image>::import( imageOrgan );
  auto imgDist = GeomHelpers::getImageDistance<CoronaryArteryTree::Image,
                                              CoronaryArteryTree::ImageDist>(img);
  double m = 0.0;
  DGtal::Z2i::Point pM;
  for(auto p: imgDist.domain()) {if (imgDist(p) > m ){m = imgDist(p); pM = p;}}
  if (verbose){
    DGtal::trace.info() << "center point found: " << pM << "maximal value:"
                        << m <<   std::endl;
  }
  constructTree(aPerf, nbTerm, imageOrgan, fgTh, verbose, pM,
                static_cast<unsigned int >(m)/2.0);
}


void
ConstructionHelpers::constructTree(double aPerf, int nbTerm,
                                   std::string imageOrgan, unsigned int fgTh,
                                   bool verbose, DGtal::Z2i::Point ptCenter,
                                   unsigned int distSearchRoot) {
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test random adds with distance constraint");
  srand (time(NULL));
  double rRoot = 1.0;
  std::string filename;
  
  CoronaryArteryTree cTree (aPerf, nbTerm, rRoot, ptCenter);
  if (imageOrgan != ""){
    auto img = DGtal::GenericReader<CoronaryArteryTree::Image>::import( imageOrgan );
    bool restrainedOK = cTree.restrainDomain(img, fgTh);
    if (restrainedOK){
      DGtal::trace.info() << "Using restrained image  " << imageOrgan << std::endl;
      cTree.myVectSegments[1].myCoordinate = ptCenter;
      cTree.myTreeCenter = ptCenter;
      if (ptCenter[0]==-1 && ptCenter[0]==-1) {
        CoronaryArteryTree::Point2D pC = cTree.getDomainCenter();
        cTree.myVectSegments[1].myCoordinate = pC;
      }
      DGtal::Z2i::Point pRoot;
      cTree.searchRootFarthest(distSearchRoot, pRoot);
      cTree.myVectSegments[0].myCoordinate = pRoot;
    }
  }
  bool isOK = false;
  unsigned int nbSeed = cTree.my_NTerm;
  for (unsigned int i = 1; i < nbSeed; i++) {
    DGtal::trace.progressBar(i, nbSeed);
    int nbSol = 0, itOpt = 0;
    CoronaryArteryTree cTreeOpt = cTree;
    double volOpt = -1.0, vol = 0.0;
    while (nbSol==0) {
      CoronaryArteryTree::Point2D pt = cTree.generateNewLocation(100);
      std::vector<unsigned int> vecN = cTree.getN_NearestSegments(pt,cTree.myNumNeighbor);
      for(size_t it=0; it<vecN.size(); it++) {
        //if(!cTree.isIntersecting(pt, cTree.findBarycenter(pt, vecN.at(it)),vecN.at(it),n))
        if(!cTree.isIntersecting(pt, cTree.findBarycenter(pt, vecN.at(it)),vecN.at(it),cTree.myNumNeighbor, 2*cTree.myVectSegments[vecN.at(it)].myRadius)) {
          CoronaryArteryTree cTree1 = cTree;
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
  cTree.myBoard.clear();
}



