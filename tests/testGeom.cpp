#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

#include "CoronaryArteryTree.h"
#include "GeomHelpers.h"
#include "ConstructionHelpers.h"

/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{
  std::string resource_dir = SAMPLE_DIR;

  ///-----------------------------------------------------------------------------------------------
  /// Test intesection on segments
  DGtal::trace.beginBlock("Testing basic intersections on two segments");
  DGtal::trace.info() << "Test intersection segments [(10, 10) (20,10)] and [(15, 20) (15,0)]";
  bool intersec1 = GeomHelpers::hasIntersection(DGtal::Z2i::RealPoint(10, 10), DGtal::Z2i::RealPoint(20, 10),
                                   DGtal::Z2i::RealPoint(15, 20), DGtal::Z2i::RealPoint(15, 0));
  DGtal::trace.info() << "Test intersection " << (intersec1 ? "OK" : "KO") << std::endl;

  DGtal::trace.info() << "Test intersection segments [(12.0, 11.0) (19,9.4)] and [(17.3, 22.0) (12,3.0)]";
  bool intersec2 = GeomHelpers::hasIntersection(DGtal::Z2i::RealPoint(12.0, 11.0), DGtal::Z2i::RealPoint(19,9.4),
                                   DGtal::Z2i::RealPoint(17.3, 22.0), DGtal::Z2i::RealPoint(12, 3.0));
  DGtal::trace.info() << "Test intersection " << (intersec2 ? "OK" : "KO") << std::endl;

  DGtal::trace.info() << "Test intersection segments [(12.0, 11.0) (19,9.4)] and [(12.0, 10.0) (19,9.4)]";
  bool intersec3 = GeomHelpers::hasIntersection(DGtal::Z2i::RealPoint(12.0, 11.0), DGtal::Z2i::RealPoint(19,9.4),
                                   DGtal::Z2i::RealPoint(12.0, 10.0), DGtal::Z2i::RealPoint(19, 9.4));
  DGtal::trace.info() << "Test intersection " << (!intersec3 ? "OK" : "KO") << std::endl;
  DGtal::trace.endBlock();
  
  
  ///-----------------------------------------------------------------------------------------------
  DGtal::trace.beginBlock("Testing basic intersections on Tree");
  CoronaryArteryTree c ( 20);
  c.addSegmentFromPoint(DGtal::Z2i::RealPoint(-10, 10), 1);
  c.addSegmentFromPoint(DGtal::Z2i::RealPoint(12, 8), 1);
  c.boardDisplay();
  c.myBoard.setPenColor(DGtal::Color::Cyan);
  c.myBoard.fillCircle(-10, 5, 1, 1);
  c.myBoard.setPenColor(DGtal::Color::Yellow);
  c.myBoard.fillCircle(-7, 6, 1, 1);
  c.myBoard.saveEPS("tmp.eps");
  c.exportBoardDisplay("testIntersection1.svg", 5.0, false);
  c.exportBoardDisplay("testIntersection1.eps", 5.0, false);
  DGtal::trace.info() << "Test intersection segments ";
  bool intersec4 = c.hasNearestIntersections(DGtal::Z2i::RealPoint(-10, -5),
                            DGtal::Z2i::RealPoint(10, -5), 10);
  DGtal::trace.info() << "Test intersection 1: no intersection " << (!intersec4 ? "OK" : "KO") << std::endl;
  bool intersec5 = c.hasNearestIntersections(DGtal::Z2i::RealPoint(-5, 5),
                            DGtal::Z2i::RealPoint(5, 5), 10);
  DGtal::trace.info() << "Test intersection 2: intersection " << (intersec5 ? "OK" : "KO") << std::endl;
  DGtal::trace.info() << "Test intersection segments ";
  bool intersec6 = c.hasNearestIntersections(c.myVectParent[1], 2,
                            DGtal::Z2i::RealPoint(-10, 6),DGtal::Z2i::RealPoint(-7, 10),  10);
  DGtal::trace.info() << "Test intersection 3:  intersection " << (!intersec6 ? "OK" : "KO") << std::endl;
  DGtal::trace.endBlock();
  ///-----------------------------------------------------------------------------------------------

  
  
  
  ///-----------------------------------------------------------------------------------------------
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test random adds with distance constraint");
  srand (time(NULL));
  CoronaryArteryTree cIntersec (DGtal::Z2i::RealPoint(0, 0), DGtal::Z2i::RealPoint(0, 30), DGtal::Z2i::RealPoint (0,  0), 100);
  for (unsigned int i = 0; i < 100; i++){
    DGtal::trace.progressBar(i, 100);
    CoronaryArteryTree::Point2D pt = cIntersec.generateNewLocation(100);
    auto nearest = cIntersec.getNearestSegment(pt);
    cIntersec.addSegmentFromPoint(pt, nearest);
  }
  cIntersec.boardDisplay();

  cIntersec.myBoard.setPenColor(DGtal::Color::Cyan);
  cIntersec.myBoard.fillCircle(0, 0, 1, 1);
  
  cIntersec.exportBoardDisplay("testNoIntersect.svg", 1.0, false);
  cIntersec.exportBoardDisplay("testNoIntersect.eps", 1.0, false);
  DGtal::trace.endBlock();

  ///-----------------------------------------------------------------------------------------------
  /// Test projection on straight line
  DGtal::trace.beginBlock("Testing projection on straight line");
  DGtal::Z2i::RealPoint p0 (0,  0);
  DGtal::Z2i::RealPoint p1 (10, 0);
  DGtal::Z2i::RealPoint p2 (5.3, 15.3);
  DGtal::Z2i::RealPoint pProj (0, 0);
  bool in = GeomHelpers::projectOnStraightLine(p0, p1, p2, pProj);
  DGtal::trace.info() << "Test intersecting segments pProj " << pProj
  << " should be (5.3, 0) " << (pProj ==  DGtal::Z2i::RealPoint(5.3, 0) && in ? "OK": "KO")  << std::endl;
  DGtal::Z2i::RealPoint pp0 (0,  0);
  DGtal::Z2i::RealPoint pp1 (3, 3);
  DGtal::Z2i::RealPoint pp2 (0, 2);
  DGtal::Z2i::RealPoint ppProj (0, 0);
  in = GeomHelpers::projectOnStraightLine(pp0, pp1, pp2, ppProj);
  DGtal::trace.info() << "Test intersecting segments ppProj " << ppProj
  << " should be (1, 1) " << (((ppProj - DGtal::Z2i::RealPoint(1.0, 1.0)).norm() < 0.000001) && in ? "OK": "KO")  << std::endl;
  DGtal::Z2i::RealPoint ppp0 (0,  0);
  DGtal::Z2i::RealPoint ppp1 (3, 3);
  DGtal::Z2i::RealPoint ppp2 (4, 4);
  DGtal::Z2i::RealPoint pppProj (0, 0);
  in = GeomHelpers::projectOnStraightLine(ppp0, ppp1, ppp2, pppProj);
  DGtal::trace.info() << "Test intersecting segments pppProj no inside " << pppProj
  << " should be in pppProj " << (((pppProj - DGtal::Z2i::RealPoint(4.0, 4.0)).norm() < 0.000001) && !in ? "OK": "KO")  << std::endl;
  DGtal::trace.endBlock();
  
    
  ///-----------------------------------------------------------------------------------------------
  /// Test projection distances
  DGtal::trace.beginBlock("Testing projection distances");
  CoronaryArteryTree ci (DGtal::Z2i::RealPoint(0, 0), DGtal::Z2i::RealPoint(0, 10), DGtal::Z2i::RealPoint (0,  0), 1);
  DGtal::Z2i::RealPoint pDirect (2,  5);
  double dist = ci.getProjDistance(1, pDirect);
  DGtal::trace.info() << "Test projection on initial segment: " << pDirect
  << " distance (should be 2) :" <<dist <<  ( dist == 2.0 ? " OK": " KO")  << std::endl;
  DGtal::Z2i::RealPoint pDirect2 (-1,  -1);
  double dist2 = ci.getProjDistance(1, pDirect2);
  DGtal::trace.info() << "Test projection on initial segment: " << pDirect
  << " distance (should be sqrt(2)) :" <<dist2 <<  ( dist2*dist2 - 2.0 < 0.0000001? " OK": " KO")  << std::endl;
  DGtal::trace.endBlock();

  
  ///-----------------------------------------------------------------------------------------------
  /// Test projection distances
  DGtal::trace.beginBlock("Testing intersections");
  std::cout << "ressource:"  << resource_dir << std::endl;
  std::stringstream ss;
  ss << resource_dir <<"shape.pgm";
  CoronaryArteryTree::Image img = DGtal::GenericReader<CoronaryArteryTree::Image>::import(ss.str());
  bool checkDomInter = GeomHelpers::checkNoIntersectDomain(img, 128, DGtal::Z2i::Point(264,196), DGtal::Z2i::Point(438,225));
  DGtal::trace.info() << "Test intersection: 264 196 and 438 225 "
  << " distance (should be true) :" << checkDomInter <<  ( checkDomInter ? " OK": " KO")  << std::endl;
  bool checkDomInter2 = GeomHelpers::checkNoIntersectDomain(img, 128, DGtal::Z2i::Point(264,196), DGtal::Z2i::Point(2,56));
  DGtal::trace.info() << "Test intersection: 264 196 and 2 56 "
  << " distance (should be false) :" << checkDomInter2 <<  ( !checkDomInter2 ? " OK": " KO")  << std::endl;
  DGtal::trace.endBlock();
  
  
  DGtal::trace.beginBlock("Testing computation of ditance map..");
  typedef typename DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> ImageDouble;
  ImageDouble imgD = GeomHelpers::getImageDistance<CoronaryArteryTree::Image, ImageDouble>(img);
  DGtal::Z2i::Point pInt(313, 201);
  DGtal::Z2i::Point pExt(20, 20);

  
  DGtal::trace.info() << " Distance map of point int "  << pInt << "distance:" <<
  (int) imgD(pInt) << " sould be > 10" << (imgD(pInt) > 10 ? " OK": "KO") << std::endl;
  
  DGtal::trace.info() << " Distance map of point int "  << pExt << "distance:" <<
  (int) imgD(pExt) << " sould be 0" << ( imgD(pExt) == 0 ? " OK": "KO") << std::endl;

  DGtal::trace.info() << "Export distance map...";
  DGtal::trace.info() << "[Done]" << std::endl;
  imgD >> "distanceMap.pgm" ;
  DGtal::trace.endBlock();

  
  DGtal::trace.beginBlock("Testing generate circle points.");
  DGtal::Z2i::DigitalSet vPt =  GeomHelpers::pointsOnCircle(DGtal::Z2i::Point(3,3), 4);
  DGtal::trace.info() << "Nb points :" << vPt.size() <<
  " Should be more than 0 "<< ( vPt.size() > 0 ? " OK": "KO") <<   std::endl;
  DGtal::trace.info() << "Should contain P(7,3) and P(-1,3) :"
  << (vPt.find(DGtal::Z2i::Point(7,3)) != vPt.end()  && vPt.find(DGtal::Z2i::Point(-1,3)) != vPt.end() ? " OK": "KO") << std::endl;
  
  
  DGtal::trace.endBlock();

  
  
  
  return EXIT_SUCCESS;
}



