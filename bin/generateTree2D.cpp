#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "CLI11.hpp"

#include "CoronaryArteryTree.h"
#include "GeomHelpers.h"
#include "ConstructionHelpers.h"




/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{
  clock_t start, end;
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  int nbTerm {500};
  double aPerf {20000};
  bool verbose {false};
  std::string nameImgDom {""}; 
  std::vector<int> postInitV {-1,-1};
  app.add_option("-n,--nbTerm,1", nbTerm, "Set the number of terminal segments.", true);
  app.add_option("-a,--aPerf,2", aPerf, "The value of the input parameter A perfusion.", true);
  app.add_option("--organDomain,-d", nameImgDom, "Define the organ domain using a mask image (organ=255).");
  auto pInit = app.add_option("-p,--posInit", postInitV, "Initial position of root, if not given the position of point is determined from the image center")
  ->expected(2);
  app.add_flag("-v,--verbose", verbose);
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------
  
  DGtal::Z2i::Point ptRoot(postInitV[0], postInitV[1]);
  CoronaryArteryTree<2> tree;
  start = clock();
  //1000 => Execution time: 129.17274900 sec
  //2000 => Execution time: 478.48590200 sec
  //3000 => Execution time: 1023.94746700 sec
  //4000 => Execution time: 1896.94450700 sec
  //5000 => Execution time: 3435.08630500 sec
  if(nameImgDom != "" && pInit->empty()){
    tree = ConstructionHelpers::constructTreeImageDomain2D<DGtal::Z2i::RealPoint>(aPerf, nbTerm, nameImgDom, 128, verbose);
  } else {
    tree = ConstructionHelpers::constructTree<2>(aPerf, nbTerm, nameImgDom, 128, verbose, ptRoot);
  }
  end = clock();
  printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);

  std::string filename;
  
  filename = "testCCO_"+std::to_string(nbTerm)+".eps";
  tree.exportBoardDisplay(filename.c_str(), 1.0);
  tree.exportBoardDisplay("result.eps", 1.0);
  tree.exportBoardDisplay("result.svg", 1.0);
  tree.myBoard.clear();
  
  return EXIT_SUCCESS;
}
