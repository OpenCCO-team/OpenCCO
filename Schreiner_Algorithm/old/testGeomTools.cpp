#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"

#include "geomhelpers.h"


/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{

  bool ok = true;
  bool inside = false;
  DGtal::trace.beginBlock("Testing project on StraightLine Test 1: (basic)");
  DGtal::Z2i::RealPoint p1 (0.0, 0.0); 
  DGtal::Z2i::RealPoint p2 (1.0, 0.0);
  DGtal::Z2i::RealPoint p3 (0.5, 1.0);
  DGtal::Z2i::RealPoint p4 (0.0, 0.0);  
  inside = projectOnStraightLine(p1, p2, p3, p4);
  DGtal::trace.info() << "p4 coordinates:"  << p4 << " Should be (0.5, 0)" <<  std::endl;
  ok = ok && ((p4 - DGtal::Z2i::RealPoint (0.5, 0.0)).norm() < 0.0000001) &&  inside ;  
  if (ok)
    DGtal::trace.info() << "TEST PASSED" << std::endl;
  else
    DGtal::trace.info() << "TEST ERRORS..." << std::endl;    
  DGtal::trace.endBlock();
  
  DGtal::trace.beginBlock("Testing project on StraightLine Test 2: (proj on diag)");
  p1 = DGtal::Z2i::RealPoint (0.0, 0.0); 
  p2 = DGtal::Z2i::RealPoint (1.0, 1.0);
  p3 = DGtal::Z2i::RealPoint (0.0, 1.0);
  p4 = DGtal::Z2i::RealPoint (0.0, 0.0);  
  inside = projectOnStraightLine(p1, p2, p3, p4);
  DGtal::trace.info() << "p4 coordinates:"  << p4 << " Should be (0.5, 0.5)" << std::endl;
  ok = ok && ((p4 - DGtal::Z2i::RealPoint (0.5, 0.5)).norm() < 0.0000001) && inside  ;
  if (ok)
    DGtal::trace.info() << "TEST PASSED" << std::endl;
  else
    DGtal::trace.info() << "TEST ERRORS..." << std::endl;    
  DGtal::trace.endBlock();

  DGtal::trace.beginBlock("Testing project on StraightLine Test 3: proj outside");
  p1 = DGtal::Z2i::RealPoint (0.0, 0.0); 
  p2 = DGtal::Z2i::RealPoint (1.0, 1.0);
  p3 = DGtal::Z2i::RealPoint (-1.0, 0.0);
  p4 = DGtal::Z2i::RealPoint (0.0, 0.0);  
  inside = projectOnStraightLine(p1, p2, p3, p4);
  DGtal::trace.info() << "p4 coordinates:"  << p4 << " Should be (-0.5, -0.5)" << std::endl;
  ok = ok && ((p4 - DGtal::Z2i::RealPoint (-0.5, -0.5)).norm() < 0.0000001) && !inside ;
  if (ok)
    DGtal::trace.info() << "TEST PASSED" << std::endl;
  else
    DGtal::trace.info() << "TEST ERRORS..." << std::endl;    
  DGtal::trace.endBlock();

  
  bool intersect = true;
  
  DGtal::trace.beginBlock("Testing intersections Test 1: Empty intersections:");
  p1 = DGtal::Z2i::RealPoint (0.0, 0.0); 
  p2 = DGtal::Z2i::RealPoint (1.0, 1.0);
  p3 = DGtal::Z2i::RealPoint (1.0, 0.5);
  p4 = DGtal::Z2i::RealPoint (1.0, -0.5);  
  DGtal::trace.info() << "Test intersection in (p1 p2):"  << p1 << " " << p2 << std::endl;
  DGtal::trace.info() << "with [p3 p4]:"  << p3 << " " << p4 << std::endl;
  intersect = hasIntersection(p1, p2, p3, p4);
  DGtal::trace.info() << "intersection ? " << (intersect? "yes": "no") << " (should be no)" << std::endl;
  if (!intersect)
    DGtal::trace.info() << "TEST PASSED" << std::endl;
  else
    DGtal::trace.info() << "TEST ERRORS..." << std::endl;    
  DGtal::trace.endBlock();

  
  DGtal::trace.beginBlock("Testing intersections Test 2: with intersections:");
  p1 = DGtal::Z2i::RealPoint (0.0, 0.0); 
  p2 = DGtal::Z2i::RealPoint (1.0, 1.0);
  p3 = DGtal::Z2i::RealPoint (0.0, 1.0);
  p4 = DGtal::Z2i::RealPoint (1.0, 0.0);  
  DGtal::trace.info() << "Test intersection in (p1 p2):"  << p1 << " " << p2 << std::endl;
  DGtal::trace.info() << "with [p3 p4]:"  << p3 << " " << p4 << std::endl;
  intersect = hasIntersection(p1, p2, p3, p4);
  DGtal::trace.info() << "intersection ? " << (intersect? "yes": "no") << " (should be yes)" << std::endl;
  if (intersect)
    DGtal::trace.info() << "TEST PASSED" << std::endl;
  else
    DGtal::trace.info() << "TEST ERRORS..." << std::endl;    
  DGtal::trace.endBlock();  


  DGtal::trace.beginBlock("Testing intersections Test 3: with parallel segments:");
  p1 = DGtal::Z2i::RealPoint (0.0, 0.0); 
  p2 = DGtal::Z2i::RealPoint (1.0, 0.0);
  p3 = DGtal::Z2i::RealPoint (0.0, 1.0);
  p4 = DGtal::Z2i::RealPoint (1.0, 1.0);  
  DGtal::trace.info() << "Test intersection in (p1 p2):"  << p1 << " " << p2 << std::endl;
  DGtal::trace.info() << "with [p3 p4]:"  << p3 << " " << p4 << std::endl;
  intersect = hasIntersection(p1, p2, p3, p4);
  DGtal::trace.info() << "intersection ? " << (intersect? "yes": "no") << " (should be no)" << std::endl;
  if (!intersect)
    DGtal::trace.info() << "TEST PASSED" << std::endl;
  else
    DGtal::trace.info() << "TEST ERRORS..." << std::endl;    
  DGtal::trace.endBlock();  



  DGtal::trace.beginBlock("Testing intersections Test 4: with coincident lines:");
  p1 = DGtal::Z2i::RealPoint (0.0, 0.0); 
  p2 = DGtal::Z2i::RealPoint (1.0, 0.0);
  p3 = DGtal::Z2i::RealPoint (2.0, 0.0);
  p4 = DGtal::Z2i::RealPoint (3.0, 0.0);  
  DGtal::trace.info() << "Test intersection in (p1 p2):"  << p1 << " " << p2 << std::endl;
  DGtal::trace.info() << "with [p3 p4]:"  << p3 << " " << p4 << std::endl;
  intersect = hasIntersection(p1, p2, p3, p4);
  DGtal::trace.info() << "intersection ? " << (intersect? "yes": "no") << " (should be no)" << std::endl;
  if (!intersect)
    DGtal::trace.info() << "TEST PASSED" << std::endl;
  else
    DGtal::trace.info() << "TEST ERRORS..." << std::endl;    
  DGtal::trace.endBlock();  

  
  
  DGtal::trace.beginBlock("Testing intersections Test 5: with various points:");
  p1 = DGtal::Z2i::RealPoint (10.0, 18.0);
  p2 = DGtal::Z2i::RealPoint (30.0, -15.3);
  p3 = DGtal::Z2i::RealPoint (8.0, -10.0);
  p4 = DGtal::Z2i::RealPoint (33.0, 022.3);
  DGtal::trace.info() << "Test intersection in (p1 p2):"  << p1 << " " << p2 << std::endl;
  DGtal::trace.info() << "with [p3 p4]:"  << p3 << " " << p4 << std::endl;
  bool intersect2 = hasIntersection(p1, p2, p3, p4);
  DGtal::trace.info() << "intersection ? " << (intersect? "yes": "no") << " (should be no)" << std::endl;
  if (intersect2)
    DGtal::trace.info() << "TEST PASSED" << std::endl;
  else
    DGtal::trace.info() << "TEST ERRORS..." << std::endl;
  DGtal::trace.endBlock();

  
  
  DGtal::trace.beginBlock("Testing intersections Test 6: test is On Right");
   p1 = DGtal::Z2i::RealPoint (0.0, 0.0);
   p2 = DGtal::Z2i::RealPoint (30.0, 20.3);
   p3 = DGtal::Z2i::RealPoint (8.0, -1.0);
   p4 = DGtal::Z2i::RealPoint (33.0, 122.3);
   DGtal::trace.info() << "Test p3 is on right (p1 p2):"  << p1 << " " << p2 << std::endl;
   DGtal::trace.info() << "with p3:"  << p3 << std::endl;
   bool isR = isOnRight(p1, p2, p3);
   DGtal::trace.info() << "isOnRight ? " << (isR? "yes": "no") << " (should be yes)" << std::endl;
   if (isR )
     DGtal::trace.info() << "TEST PASSED" << std::endl;
   else
     DGtal::trace.info() << "TEST ERRORS..." << std::endl;
   DGtal::trace.endBlock();

  return EXIT_SUCCESS;
}
