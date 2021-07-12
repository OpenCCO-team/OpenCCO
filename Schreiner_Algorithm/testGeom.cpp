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
  
  return EXIT_SUCCESS;
}



