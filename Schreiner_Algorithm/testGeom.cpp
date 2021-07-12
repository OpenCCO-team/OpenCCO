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
  
  DGtal::trace.beginBlock("Testing basic intersections on two segments");
  DGtal::trace.info() << "Test intersection segments [(10, 10) (20,10)] and [(15, 20) (15,0)]";
  bool intersec1 = hasIntersection(DGtal::Z2i::RealPoint(10, 10), DGtal::Z2i::RealPoint(20, 10),
                                   DGtal::Z2i::RealPoint(15, 20), DGtal::Z2i::RealPoint(15, 0));
  DGtal::trace.info() << "Test intersection " << (intersec1 ? "OK" : "KO") << std::endl;

  DGtal::trace.info() << "Test intersection segments [(12.0, 11.0) (19,9.4)] and [(17.3, 22.0) (12,3.0)]";
  bool intersec2 = hasIntersection(DGtal::Z2i::RealPoint(12.0, 11.0), DGtal::Z2i::RealPoint(19,9.4),
                                   DGtal::Z2i::RealPoint(17.3, 22.0), DGtal::Z2i::RealPoint(12, 3.0));
  DGtal::trace.info() << "Test intersection " << (intersec2 ? "OK" : "KO") << std::endl;

  DGtal::trace.info() << "Test intersection segments [(12.0, 11.0) (19,9.4)] and [(12.0, 10.0) (19,9.4)]";
  bool intersec3 = hasIntersection(DGtal::Z2i::RealPoint(12.0, 11.0), DGtal::Z2i::RealPoint(19,9.4),
                                   DGtal::Z2i::RealPoint(12.0, 10.0), DGtal::Z2i::RealPoint(19, 9,4));
  DGtal::trace.info() << "Test intersection " << (!intersec3 ? "OK" : "KO") << std::endl;
  DGtal::trace.endBlock();
  
  
  DGtal::trace.beginBlock("Testing basic intersections on Tree");
  CoronaryArteryTree c (DGtal::Z2i::RealPoint(0, 0), 2000, 1);
  c.addSegmentFromPoint(DGtal::Z2i::RealPoint(-10, 10), 1);
  c.addSegmentFromPoint(DGtal::Z2i::RealPoint(12, 8), 1);
  c.boardDisplay();
  c.myBoard.setFillColor(DGtal::Color::Cyan);
  c.myBoard.drawCircle(-10, 5, 1, 1);
  c.myBoard.setFillColor(DGtal::Color::Yellow);
  c.myBoard.drawCircle(-7, 6, 1, 1);
  c.exportBoardDisplay("testIntersection1.svg", false);
  c.exportBoardDisplay("testIntersection1.eps", false);
  DGtal::trace.info() << "Test intersection segments ";
  bool intersec4 = c.hasNearestIntersections(DGtal::Z2i::RealPoint(-10, -5),
                            DGtal::Z2i::RealPoint(10, -5), 10);
  DGtal::trace.info() << "Test intersection 1: no intersection " << (!intersec4 ? "OK" : "KO") << std::endl;
  bool intersec5 = c.hasNearestIntersections(DGtal::Z2i::RealPoint(-5, 5),
                            DGtal::Z2i::RealPoint(5, 5), 10);
  DGtal::trace.info() << "Test intersection 2: intersection " << (intersec5 ? "OK" : "KO") << std::endl;
  DGtal::trace.info() << "Test intersection segments ";
  bool intersec6 = c.hasNearestIntersections(c.myVectParent[2], 2,
                            DGtal::Z2i::RealPoint(-10, 6),DGtal::Z2i::RealPoint(-7, 10),  10);
  DGtal::trace.info() << "Test intersection 3:  intersection " << (!intersec6 ? "OK" : "KO") << std::endl;
  DGtal::trace.endBlock();

  
  
  
  
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test random adds with distance constraint");
  srand (time(NULL));
  CoronaryArteryTree cIntersec (DGtal::Z2i::RealPoint(0, 250), 20000000, 100);
  
  for (unsigned int i = 0; i < 500; i++){
    DGtal::trace.progressBar(i, 500);
    CoronaryArteryTree::Point2D pt = cIntersec.generateNewLocation(100);
    auto nearest = cIntersec.getNearestSegment(pt);
    cIntersec.addSegmentFromPoint(pt, nearest, 1.0, 1.0);
    cIntersec.udpatePerfusionArea();
    cIntersec.updateTreshold();
  }
  cIntersec.exportBoardDisplay("testNoIntersect.svg", true);
  cIntersec.exportBoardDisplay("testNoIntersect.eps", true);
  
  
  
  return EXIT_SUCCESS;
}



