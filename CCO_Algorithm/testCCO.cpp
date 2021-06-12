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
  
  unsigned int nbSeed = 50;
  for (unsigned int i = 0; i < nbSeed; i++){
    DGtal::trace.progressBar(i, nbSeed);
    CoronaryArteryTree::Point2D pt = cTree.generateNewLocation(100);
    auto nearest = cTree.getNearestSegment(pt);
    cTree.addSegmentFromPoint(pt, nearest, 1.0, 1.0);
    cTree.udpatePerfusionArea();
    cTree.updateTreshold();
  }
  
  cTree.exportBoardDisplay("testCCO1.eps", true);
  cTree.myBoard.clear();
  
  unsigned int idSeg = nbSeed - 1;
  unsigned int idTarget = cTree.myVectParent[idSeg];
  
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
