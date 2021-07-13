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
  ///-----------------------------------------------------------------------------------------------
  /// Test intesection on segments
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
                                   DGtal::Z2i::RealPoint(12.0, 10.0), DGtal::Z2i::RealPoint(19, 9.4));
  DGtal::trace.info() << "Test intersection " << (!intersec3 ? "OK" : "KO") << std::endl;
  DGtal::trace.endBlock();
  
  
  ///-----------------------------------------------------------------------------------------------
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
  ///-----------------------------------------------------------------------------------------------

  
  
  
  ///-----------------------------------------------------------------------------------------------
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
  cIntersec.boardDisplay();

  cIntersec.myBoard.setFillColor(DGtal::Color::Cyan);
  cIntersec.myBoard.drawCircle(0, 0, 50, 1);

  cIntersec.exportBoardDisplay("testNoIntersect.svg", false);
  cIntersec.exportBoardDisplay("testNoIntersect.eps", false);
  DGtal::trace.endBlock();

  ///-----------------------------------------------------------------------------------------------
  /// Test projection on straight line
  DGtal::trace.beginBlock("Testing projection on straight line");
  DGtal::Z2i::RealPoint p0 (0,  0);
  DGtal::Z2i::RealPoint p1 (10, 0);
  DGtal::Z2i::RealPoint p2 (5.3, 15.3);
  DGtal::Z2i::RealPoint pProj (0, 0);
  bool in = projectOnStraightLine(p0, p1, p2, pProj);
  DGtal::trace.info() << "Test intersecting segments pProj " << pProj
  << " should be (5.3, 0) " << (pProj ==  DGtal::Z2i::RealPoint(5.3, 0) && in ? "OK": "KO")  << std::endl;
  DGtal::Z2i::RealPoint pp0 (0,  0);
  DGtal::Z2i::RealPoint pp1 (3, 3);
  DGtal::Z2i::RealPoint pp2 (0, 2);
  DGtal::Z2i::RealPoint ppProj (0, 0);
  in = projectOnStraightLine(pp0, pp1, pp2, ppProj);
  DGtal::trace.info() << "Test intersecting segments ppProj " << ppProj
  << " should be (1, 1) " << (((ppProj - DGtal::Z2i::RealPoint(1.0, 1.0)).norm() < 0.000001) && in ? "OK": "KO")  << std::endl;
  DGtal::Z2i::RealPoint ppp0 (0,  0);
  DGtal::Z2i::RealPoint ppp1 (3, 3);
  DGtal::Z2i::RealPoint ppp2 (4, 4);
  DGtal::Z2i::RealPoint pppProj (0, 0);
  in = projectOnStraightLine(ppp0, ppp1, ppp2, pppProj);
  DGtal::trace.info() << "Test intersecting segments pppProj no inside " << pppProj
  << " should be in pppProj " << (((pppProj - DGtal::Z2i::RealPoint(4.0, 4.0)).norm() < 0.000001) && !in ? "OK": "KO")  << std::endl;
  DGtal::trace.endBlock();
  
    
  ///-----------------------------------------------------------------------------------------------
  /// Test projection distances
  DGtal::trace.beginBlock("Testing projection distances");
  CoronaryArteryTree ci (DGtal::Z2i::RealPoint(0, 0), 2000, 1);
  
  DGtal::Z2i::RealPoint pDirect (2,  5);
  double dist = ci.getProjDistance(1, pDirect);
  DGtal::trace.info() << "Test projection on initial segment: " << pDirect
  << " distance (should be 2) :" <<dist <<  ( dist == 2.0 ? " OK": " KO")  << std::endl;
  DGtal::Z2i::RealPoint pDirect2 (-1,  -1);
  double dist2 = ci.getProjDistance(1, pDirect2);
  DGtal::trace.info() << "Test projection on initial segment: " << pDirect
  << " distance (should be sqrt(2)) :" <<dist2 <<  ( dist2*dist2 - 2.0 < 0.0000001? " OK": " KO")  << std::endl;
  DGtal::trace.endBlock();

  
  return EXIT_SUCCESS;
}



