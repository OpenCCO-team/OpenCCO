#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "CLI11.hpp"

#include "CoronaryArteryTree.h"
#include "GeomHelpers.h"
#include "ExpandTreeHelpers.h"



/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{
  std::string resource_dir = SAMPLE_DIR;
  clock_t start, end;
  
  int nbTerm {300};
  double aPerf {20000};
  bool verbose {false};
  std::stringstream ss;
  ss << resource_dir <<"shape.png";
  std::string nameImgDom  = ss.str();
  //----------------------
  // Example to construct a simple tree in 2D from mask respresented by the file (shape.png)
  //----------------------
  start = clock();
  // 1. Type definition, of domain controller and tree.
  typedef ImageMaskDomainCtrl<2> TImgContrl;
  typedef  CoronaryArteryTree<TImgContrl, 2> TTreeMaskDom;
  
  // 2. Domain controller construction
  TImgContrl aDomCtr = TImgContrl(nameImgDom, 128, 100);
  
  // 3. Tree construction using center
  TTreeMaskDom tree  (aPerf, nbTerm, aDomCtr);

  // 5. Tree construction
  ExpandTreeHelpers::initFirstElemTree(tree, true);
  ExpandTreeHelpers::expandTree(tree, true);
  
  // 5. exporting the result
  tree.exportBoardDisplay("testResult2TreeMaskedDomain.eps", 1.0);
  end = clock();
  printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
  
  //----------------------
  // Example to construct a simple tree in 2D from implicit domain
  //----------------------

  // Example to construct a simple tree in 2D from mask respresented by the file (shape.png)
  // 1. Type definition, of domain controller and tree.
  typedef CircularDomainCtrl<2> ImplicitContrl;
  typedef CoronaryArteryTree<ImplicitContrl, 2> TTreeCircDom;

  // 2. Domain controller construction
  ImplicitContrl aDomCtrImp(1.0, TImgContrl::TPoint(0,0));

  // 3. Tree construction using center
  TTreeCircDom treeImpl (aPerf, nbTerm, aDomCtrImp);

  // 4. Tree expansion
  ExpandTreeHelpers::initFirstElemTree(treeImpl, true);
  ExpandTreeHelpers::expandTree(treeImpl, true);

  // 5. Tree export
  treeImpl.exportBoardDisplay("testResult2TreeMaskedDomain.eps", 1.0);


  return EXIT_SUCCESS;
}
