#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "CLI11.hpp"

#include "CoronaryArteryTree.h"
#include "geomhelpers.h"
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
   int nbTerm {300};
   double aPerf {20000};
   app.add_option("-n,--nbTerm,1", nbTerm, "Set the number of terminal segments.", true);
   app.add_option("-a,--aPerf,1", aPerf, "The value of the input parameter A perfusion.", true);

   
   start = clock();
  //1000 => Execution time: 129.17274900 sec
  //2000 => Execution time: 478.48590200 sec
  //3000 => Execution time: 1023.94746700 sec
  ConstructionHelpers::constructTree(aPerf, nbTerm);
  end = clock();
  printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
  
  return EXIT_SUCCESS;
}
