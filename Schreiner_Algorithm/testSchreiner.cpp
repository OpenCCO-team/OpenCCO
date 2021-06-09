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

  bool ok = true;
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: (base)");
  CoronaryArteryTree c (DGtal::Z2i::RealPoint(0, 25), 2000, 1);
  DGtal::trace.info() << c;
  c.exportDisplay("testBase0.svg");
  c.exportDisplay("testBase0.eps");
  DGtal::trace.endBlock();

//   DGtal::trace.beginBlock("Testing class CoronaryArteryTree: construction base (1)");
//   c.addFirstSegment(DGtal::Z2i::RealPoint(10, 10));
//   DGtal::Z2i::RealPoint p = c.getSegmentCenter(1);
//   DGtal::trace.info() <<"Center point of first segment: " <<  p << "(should be (5, 17)" << std::endl;
//   ok = ok && p == DGtal::Z2i::RealPoint(5, 17);
//   if (ok)
//     DGtal::trace.info() << "TEST PASSED" << std::endl;
//   else
//     DGtal::trace.info() << "TEST ERRORS..." << std::endl;
//   DGtal::trace.endBlock();

  
  
  
//   DGtal::trace.beginBlock("Testing class CoronaryArteryTree: construction base (2)");
//   bool ok2 = false;
//   auto nearest = c.getNearestSegment(DGtal::Z2i::RealPoint(15,15));
//   DGtal::trace.info() <<"Nearest Segment of (15, 15): " <<  c.myVectSegments[nearest].myCoordinate << " (should be (10, 10)" << std::endl;
//   ok2 = c.myVectSegments[nearest].myCoordinate == DGtal::Z2i::RealPoint(10,10);
//   if (ok2)
//      DGtal::trace.info() << "TEST PASSED" << std::endl;
//    else
//      DGtal::trace.info() << "TEST ERRORS..." << std::endl;
//   DGtal::trace.endBlock();

  
  
//   DGtal::trace.beginBlock("Testing class CoronaryArteryTree: construction base (3)");
//   bool ok3 = false;
//   c.addSegmentFromPoint(DGtal::Z2i::RealPoint(15, 20));
//   unsigned int addIndex = c.myVectSegments.size()-1;
//   DGtal::trace.info() <<"Parent of new segment " <<  c.myVectSegments[c.myVectParent[addIndex]].myCoordinate << " (should be (5, 17.5)" << std::endl;
//   ok3 = c.myVectSegments[c.myVectParent[addIndex]].myCoordinate == DGtal::Z2i::RealPoint(5,17.5);
//   unsigned int indexFirstDaugther = c.myVectDaughters[c.myVectParent[addIndex]].first;
//   unsigned int indexSecDaugther = c.myVectDaughters[c.myVectParent[addIndex]].second;

//   DGtal::trace.info() <<"First Daugther of parent of new segment " <<  c.myVectSegments[indexFirstDaugther].myCoordinate  << " (should be (10, 10)" << std::endl;
//   ok3 = ok3 &&  c.myVectSegments[indexFirstDaugther].myCoordinate == DGtal::Z2i::RealPoint(10, 10);

//   DGtal::trace.info() <<"Second Daugther of parent of new segment " <<  c.myVectSegments[indexSecDaugther].myCoordinate  << " (should be (15, 20)" << std::endl;
//   ok3 = ok3 &&  c.myVectSegments[indexSecDaugther].myCoordinate == DGtal::Z2i::RealPoint(15, 20);

  
//   if (ok3)
//       DGtal::trace.info() << "TEST PASSED" << std::endl;
//     else
//       DGtal::trace.info() << "TEST ERRORS..." << std::endl;
//    DGtal::trace.endBlock();

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
