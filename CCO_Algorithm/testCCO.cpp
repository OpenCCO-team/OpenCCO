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
  CoronaryArteryTree cTree (DGtal::Z2i::RealPoint(0, 250), 20000000, 100);
  
  for (unsigned int i = 0; i < 50; i++){
    DGtal::trace.progressBar(i, 50);
    CoronaryArteryTree::Point2D pt = cTree.generateNewLocation(100);
    cTree.udpatePerfusionArea();
    cTree.updateTreshold();
  }
  
  cTree.exportBoardDisplay("testCCO.svg", true);
  cTree.exportBoardDisplay("testCCO.eps", true);
  
  DGtal::trace.endBlock();
 
  
  return EXIT_SUCCESS;
}
