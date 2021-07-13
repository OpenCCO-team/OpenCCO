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
  CoronaryArteryTree cTree (pRoot, 2000000, 10);
  
  std::string filename;
  unsigned int nbSeed = cTree.my_NTerm;
  bool isOK;
  for (unsigned int i = 0; i < nbSeed; i++) {
    //DGtal::trace.progressBar(i, nbSeed);
    isOK = false;
    while (!isOK) {
      CoronaryArteryTree::Point2D pt = cTree.generateNewLocation(100);
      auto nearest = cTree.getNearestSegment(pt);
      /*
      filename = "testCCO_"+std::to_string(nearest)+"A.eps";
      cTree.exportBoardDisplay(filename.c_str(), true);
      cTree.myBoard.clear();
      */
      isOK = cTree.isAddable(pt,nearest, 100);
      /*
      std::cout<<"isOK="<<isOK<<std::endl;
      if(isOK) {
        filename = "testCCO_"+std::to_string(nearest)+"D.eps";
        cTree.exportBoardDisplay(filename.c_str(), true);
        cTree.myBoard.clear();
      }
      */
      //std::cout<<"isOK="<<isOK<<std::endl;
      if(isOK) {
        filename = "testCCO_"+std::to_string(i)+".eps";
        cTree.exportBoardDisplay(filename.c_str(), true);
        cTree.myBoard.clear();
      }
    }
    //std::cout<<i<<" => Total volume="<<cTree.computeTotalVolume(1)<<std::endl;
  }
    
  unsigned int n = 5, nbSol = 0, itOpt = 0;
  CoronaryArteryTree cTreeOpt = cTree;
  double volOpt = -1.0, vol = 0.0;
  for (unsigned int i = 0; i < 1; i++){
    //DGtal::trace.progressBar(i, nbSeed);
    while (nbSol==0) {
      CoronaryArteryTree::Point2D pt = cTree.generateNewLocation(100);
      std::vector<unsigned int> vecN = cTree.getN_NearestSegments(pt,n);
      for(size_t it=0; it<vecN.size(); it++) {
        CoronaryArteryTree cTree1 = cTree;
        isOK = cTree1.isAddable(pt,vecN.at(it), 100);
        if(isOK) {
          vol = cTree1.computeTotalVolume(1);
          std::cout<<" => Test Total volume ("<<i<<"-"<<it<<")="<<vol<<std::endl;
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
          filename = "testCCO_V"+std::to_string(i)+"_"+std::to_string(it)+".eps";
          cTree1.exportBoardDisplay(filename.c_str(), true);
          cTree1.myBoard.clear();
          nbSol++;
        }
      }
    }
  }
  
  std::cout<<"itOpt="<<itOpt<<std::endl;
  
  cTree.exportBoardDisplay("testCCO1.eps", true);
  cTree.myBoard.clear();
  
  cTreeOpt.exportBoardDisplay("testCCO2.eps", true);
  cTreeOpt.myBoard.clear();
  
  return EXIT_SUCCESS;
}
