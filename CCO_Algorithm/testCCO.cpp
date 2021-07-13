#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"


#include "CoronaryArteryTree.h"
#include "geomhelpers.h"

/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test random adds with distance constraint");
  srand (time(NULL));
  DGtal::Z2i::RealPoint pRoot;
  double aPerf = 2000000;
  int nbTerm = 100;
  double rRoot = 1.0;//1.0/nbTerm;
  
  CoronaryArteryTree cTree (pRoot, aPerf, nbTerm, rRoot);
  
  unsigned int n = 10;
  bool isOK = false;
  std::string filename;
  unsigned int nbSeed = cTree.my_NTerm - 1;
  for (unsigned int i = 0; i < nbSeed; i++) {
    //DGtal::trace.progressBar(i, nbSeed);
    int nbSol = 0, itOpt = 0;
    CoronaryArteryTree cTreeOpt = cTree;
    double volOpt = -1.0, vol = 0.0;
    while (nbSol==0) {
      CoronaryArteryTree::Point2D pt = cTree.generateNewLocation(100);
      std::vector<unsigned int> vecN = cTree.getN_NearestSegments(pt,n);
      for(size_t it=0; it<vecN.size(); it++) {
        if(!cTree.isIntersecting(pt,vecN.at(it),n)) {
          CoronaryArteryTree cTree1 = cTree;
          isOK = cTree1.isAddable(pt,vecN.at(it), 100);
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
    /*
    filename = "testCCO_V_A"+std::to_string(i)+".eps";
    cTreeOpt.exportBoardDisplay(filename.c_str(), true);
    cTreeOpt.myBoard.clear();
    */
    cTree = cTreeOpt;
    //cTree.updateScale(sqrt(1.0+(1.0/(i+1.0))));
    std::cout<<"it="<<i<<"=> Aperf="<<cTree.myRsupp*cTree.myRsupp*M_PI<<std::endl;
    /*
    filename = "testCCO_V_B"+std::to_string(i)+".eps";
    cTreeOpt.exportBoardDisplay(filename.c_str(), true);
    cTreeOpt.myBoard.clear();
    */
    
  }
  std::cout<<"====> Aperf="<<cTree.myRsupp*cTree.myRsupp*M_PI<<" == "<<aPerf<<std::endl;

  cTree.exportBoardDisplay("testCCO1.eps", true);
  cTree.myBoard.clear();
  
  return EXIT_SUCCESS;
}
