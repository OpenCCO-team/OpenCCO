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
  CoronaryArteryTree cTree (DGtal::Z2i::RealPoint(0, 250), 200000, 100);
  
  for (unsigned int i = 0; i < 2; i++){
    DGtal::trace.progressBar(i, 2);
    CoronaryArteryTree::Point2D pt = cTree.generateNewLocation(100);
    auto nearest = cTree.getNearestSegment(pt);
    cTree.addSegmentFromPoint(pt, nearest, 1.0, 1.0);
    cTree.udpatePerfusionArea();
    cTree.updateTreshold();
  }
  
 
  
  cTree.myBoard.setPenColor(DGtal::Color::Yellow);
  cTree.myBoard.setFillColor(DGtal::Color::Yellow);
  cTree.myBoard.setLineWidth(1);
  unsigned int idTarget = cTree.myVectParent[2];
  cTree.myBoard.drawCircle(cTree.myVectSegments[idTarget].myCoordinate[0],
                           cTree.myVectSegments[idTarget].myCoordinate[1], 1, 0);
  cTree.myBoard.setLineWidth(50);

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

  cTree.myBoard.setFillColor(DGtal::Color::None);

  //cTree.exportBoardDisplay("testCCO1.svg", true);
  cTree.exportBoardDisplay("testCCO1.eps", true);
  
  cTree.kamyiaOptimization(idTarget);
  cTree.myBoard.setFillColor(DGtal::Color::None);

  //cTree.exportBoardDisplay("testCCO2.svg", true);
  cTree.exportBoardDisplay("testCCO2.eps", true);

  DGtal::trace.endBlock();
 
  
  return EXIT_SUCCESS;
}
