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
  
  DGtal::trace.beginBlock("esting base constructor with initial segment (should give a first random segment).");
  CoronaryArteryTree c (DGtal::Z2i::RealPoint(0, 25), 2000, 1);
  DGtal::trace.info() << c;
  c.exportBoardDisplay("testBase0.svg", true);
  c.exportBoardDisplay("testBase0.eps", true);
  DGtal::trace.endBlock();
  
  
  
  
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: construction base (1)!!!");
  bool ok = true;
  auto nearest = c.getNearestSegment(DGtal::Z2i::RealPoint(15,15));
  DGtal::trace.info() <<"Nearest Segment of (15, 15): " <<  c.myVectSegments[nearest].myCoordinate << " (should be " <<  c.myVectSegments[1].myCoordinate  << std::endl;
  ok = c.myVectSegments[nearest].myCoordinate == c.myVectSegments[1].myCoordinate;
  if (ok)
    DGtal::trace.info() << "TEST PASSED" << std::endl;
  else
    DGtal::trace.info() << "TEST ERRORS..." << std::endl;
  DGtal::trace.endBlock();
  
  
  
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: construction base (2)");
  bool ok2 = false;
  auto pRan = c.myVectSegments[1].myCoordinate;
  c.addSegmentFromPoint(DGtal::Z2i::RealPoint(15, 20), c.myVectSegments[1].myRadius);
  c.exportBoardDisplay("testBase2.svg", true);
  c.exportBoardDisplay("testBase2.eps", true);
  
  unsigned int addIndex = c.myVectSegments.size()-1;
  DGtal::trace.info() <<"Parent of new segment "
  << c.myVectSegments[c.myVectParent[addIndex]].myCoordinate << " (should be "
  << c.myVectSegments[1].myCoordinate << ")" << std::endl;
  ok2 = c.myVectSegments[c.myVectParent[addIndex]].myCoordinate == c.myVectSegments[1].myCoordinate;
  
  unsigned int indexLeft = c.myVectChildren[c.myVectParent[addIndex]].first;
  unsigned int indexRight = c.myVectChildren[c.myVectParent[addIndex]].second;
  
  DGtal::trace.info() <<"First Children of parent of new segment " <<  c.myVectSegments[indexLeft].myCoordinate  << " (should be " << pRan << std::endl;
  ok2 = ok2 &&  c.myVectSegments[indexLeft].myCoordinate == pRan;
  
  DGtal::trace.info() <<"Second Children of parent of new segment " <<  c.myVectSegments[indexRight].myCoordinate  << " (should be (15, 20)" << std::endl;
  ok2 = ok2 &&  c.myVectSegments[indexRight].myCoordinate == DGtal::Z2i::RealPoint(15, 20);
  
  
  if (ok2)
    DGtal::trace.info() << "TEST PASSED" << std::endl;
  else
    DGtal::trace.info() << "TEST ERRORS..." << std::endl;
  DGtal::trace.endBlock();
  
  
  
  
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: display constructions steps (3)");
  nearest = c.getNearestSegment(DGtal::Z2i::RealPoint(-2,6));
  c.addSegmentFromPoint(DGtal::Z2i::RealPoint(-2, 6), nearest, 1.0, 1.0);
  c.exportBoardDisplay("testBase3.svg", true);
  c.exportBoardDisplay("testBase3.eps", true);
  nearest = c.getNearestSegment(DGtal::Z2i::RealPoint(15,-5));
  c.addSegmentFromPoint(DGtal::Z2i::RealPoint(15, -5), nearest, 0.5, 0.3);
  c.exportBoardDisplay("testBase4.svg", true);
  c.exportBoardDisplay("testBase4.eps", true);
  DGtal::trace.endBlock();
  
  
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test random adds");
  srand (time(NULL));
  CoronaryArteryTree cRand (DGtal::Z2i::RealPoint(0, 250), 2000, 1);
  
  for (unsigned int i = 0; i < 1000; i++){
    CoronaryArteryTree::Point2D pt = generateRandomPtOnDisk(cRand.myTreeCenter, cRand.my_rPerf);
    
    nearest = cRand.getNearestSegment(pt);
    cRand.addSegmentFromPoint(pt, nearest, 1.0, 1.0);
  }
  
  cRand.exportBoardDisplay("testRandomAdd.svg", true);
  cRand.exportBoardDisplay("testRandomAdd.eps", true);
  
  cRand.boardDisplay();
  cRand.exportBoardDisplay("toto.eps", false);
  DGtal::trace.endBlock();
  
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test random adds");
  std::vector<unsigned int> v = c.getPathToRoot(c.myVectSegments[c.myVectSegments.size()-2]);
  c.myBoard.setLineWidth(20.0);
  c.myBoard.setPenColor(DGtal::Color::Yellow);
  c.myBoard.setFillColor(DGtal::Color::Yellow);

  for (unsigned int u : v){
    auto seg = c.myVectSegments[u];
    auto seg2 = c.myVectSegments[c.myVectParent[seg.myIndex]];
    c.myBoard.drawLine(seg.myCoordinate[0], seg.myCoordinate[1],
                       seg2.myCoordinate[0], seg2.myCoordinate[1], 0);
  }
  c.exportBoardDisplay("testBase5.svg", false);
  c.exportBoardDisplay("testBase5.eps", false);

  DGtal::trace.endBlock();

  
  
  return EXIT_SUCCESS;
}
