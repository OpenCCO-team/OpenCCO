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
  c.exportDisplay("testBase0.svg");
  c.exportDisplay("testBase0.eps");
  DGtal::trace.endBlock();
  
  
  
  
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: construction base (1)");
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
  c.addSegmentFromPoint(DGtal::Z2i::RealPoint(15, 20));
  c.exportDisplay("testBase2.svg");
  c.exportDisplay("testBase2.eps");
  
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
  
  //   c.addSegmentFromPoint(DGtal::Z2i::RealPoint(5, 3));
  //   c.addSegmentFromPoint(DGtal::Z2i::RealPoint(9, -2));
  
  
  //   DGtal::trace.beginBlock("Testing class CoronaryArteryTree: display constructions steps");
  //   c.exportDisplay("testStep1.svg");
  //   c.exportDisplay("testStep1.eps");
  //   DGtal::trace.endBlock();
  
  //   DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test random adds");
  //   srand (time(NULL));
  //   CoronaryArteryTree cRand (DGtal::Z2i::RealPoint(0, 25), 2000, 1);
  
  //   cRand.addFirstSegment(DGtal::Z2i::RealPoint(0, 0));
  
  //   for (unsigned int i = 0; i < 100; i++){
  //     double x = rand() % ((int)(c.myRsupp*200.0));
  //     double y = rand() % ((int)(c.myRsupp*200.0));
  //     x = (x/100.0 - c.myRsupp);
  //     y = (y/100.0 - c.myRsupp);
  //     if (isInsideCircle(DGtal::Z2i::RealPoint(0, 0), DGtal::Z2i::RealPoint(x, y),  c.myRsupp)){
  
  //       cRand.addSegmentFromPointWithBarycenter(DGtal::Z2i::RealPoint(x, y));
  //     }
  //   }
  
  
  // //    cout << " le point aleatoire est simule. Ses coordonnées sont: x = " << x << " et y = " << y <<endl;
  
  
  //   cRand.exportDisplay("testRandomAdd.svg");
  //   cRand.exportDisplay("testRandomAdd.eps");
  //   DGtal::trace.endBlock();
  
  
  return EXIT_SUCCESS;
}
