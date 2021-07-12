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
  CoronaryArteryTree cTree (pRoot, 2000000, 100);
  
  std::string filename;
  unsigned int nbSeed = 100;
  bool isOK;
  for (unsigned int i = 0; i < nbSeed; i++){
    DGtal::trace.progressBar(i, nbSeed);
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
    }
    //cTree.addSegmentFromPoint(pt, nearest, 1.0, 1.0);
    /*
    filename = "testCCO_"+std::to_string(i+1)+"A.eps";
    cTree.exportBoardDisplay(filename.c_str(), true);
    cTree.myBoard.clear();
    cTree.kamyiaOptimization(nearest);
    filename = "testCCO_"+std::to_string(i+1)+"B.eps";
    cTree.exportBoardDisplay(filename.c_str(), true);
    cTree.myBoard.clear();
    //cTree.udpatePerfusionArea();
    //cTree.updateTreshold();
    */
  }
  
  
  cTree.exportBoardDisplay("testCCO1.eps", true);
  cTree.myBoard.clear();
  
  unsigned int idSeg = nbSeed;//nbSeed - 1;
  unsigned int idTarget = idSeg;//cTree.myVectParent[idSeg];
  
  cTree.myBoard.setPenColor(DGtal::Color::Green);
  //cTree.myBoard.setFillColor(DGtal::Color::Yellow);
  cTree.myBoard.setLineWidth(1);
  
  cTree.myBoard.drawCircle(cTree.myVectSegments[idTarget].myCoordinate[0],
                           cTree.myVectSegments[idTarget].myCoordinate[1], 1, 0);
  //cTree.myBoard.setLineWidth(20);

  cTree.myBoard.drawLine(cTree.myVectSegments[idTarget].myCoordinate[0],
                           cTree.myVectSegments[idTarget].myCoordinate[1],
                         cTree.myVectSegments[cTree.myVectParent[idTarget]].myCoordinate[0],
                          cTree.myVectSegments[cTree.myVectParent[idTarget]].myCoordinate[1],
                         0);
  
  cTree.myBoard.drawLine(cTree.myVectSegments[idTarget].myCoordinate[0],
                           cTree.myVectSegments[idTarget].myCoordinate[1],
                         cTree.myVectSegments[cTree.myVectChildren[idTarget].first].myCoordinate[0],
                          cTree.myVectSegments[cTree.myVectChildren[idTarget].first].myCoordinate[1],
                         0);
  cTree.myBoard.drawLine(cTree.myVectSegments[idTarget].myCoordinate[0],
                           cTree.myVectSegments[idTarget].myCoordinate[1],
                         cTree.myVectSegments[cTree.myVectChildren[idTarget].second].myCoordinate[0],
                          cTree.myVectSegments[cTree.myVectChildren[idTarget].second].myCoordinate[1],
                         0);
  cTree.myBoard.setLineWidth(1);

  cTree.kamyiaOptimization(idTarget);
  
  //cTree.exportBoardDisplay("testCCO2.svg", true);
  cTree.exportBoardDisplay("testCCO2.eps", true);

  DGtal::trace.endBlock();
 
  
  return EXIT_SUCCESS;
}
