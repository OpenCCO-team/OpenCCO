
#include "ConstructionHelpers.h"
#include "CoronaryArteryTree.h"





void
ConstructionHelpers::constructTree(double aPerf, int nbTerm) {
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test random adds with distance constraint");
  srand (time(NULL));
  double rRoot = 1.0;//10.0/nbTerm;
  std::string filename;
  
  CoronaryArteryTree cTree (aPerf, nbTerm, rRoot);
  
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
          isOK = cTree1.isAddable(pt,vecN.at(it), 100, 0.01, cTree1.myNumNeighbor);
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
  std::cout<<"====> Aperf="<<cTree.myRsupp*cTree.myRsupp*cTree.my_NTerm*M_PI<<" == "<<aPerf<<std::endl;
  
  filename = "testCCO_"+std::to_string(nbTerm)+".eps";
  cTree.exportBoardDisplay(filename.c_str(), 1.0);
  cTree.myBoard.clear();
}
