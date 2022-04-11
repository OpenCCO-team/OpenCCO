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

#include "DGtal/io/viewers/Viewer3D.h"

/**
 * @brief main function call
 *
 */
int main(int argc, char **argv)
{
  QApplication application(argc,argv);
  
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
  
  DGtal::Z3i::Point ptRoot(postInitV[0], postInitV[1], 0);
  CoronaryArteryTree<3> tree;
  start = clock();
  //1000 => Execution time: 129.17274900 sec
  //2000 => Execution time: 478.48590200 sec
  //3000 => Execution time: 1023.94746700 sec
  //4000 => Execution time: 1896.94450700 sec
  //5000 => Execution time: 3435.08630500 sec
  if(nameImgDom != "" && pInit->empty()){
    tree = ConstructionHelpers::constructTreeImageDomain3D<DGtal::Z3i::RealPoint>(aPerf, nbTerm, nameImgDom, 128, verbose);
  } else {
    tree = ConstructionHelpers::constructTree<3>(aPerf, nbTerm, nameImgDom, 128, verbose, ptRoot);
  }
  end = clock();
  printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);

  typedef DGtal::Viewer3D<> MyViewer;
  MyViewer viewer;
  viewer.show();
  unsigned int i = 0;
  double thickness = 1;
  viewer << DGtal::CustomColors3D(DGtal::Color(0,0,250),DGtal::Color(0,0,250));
  DGtal::Z3i::RealPoint p1 = tree.myVectSegments[1].myCoordinate;
  DGtal::Z3i::RealPoint p2 = tree.myVectSegments[tree.myVectParent[tree.myVectSegments[1].myIndex]].myCoordinate;
  viewer.addBall(p2,tree.myVectSegments[1].myRadius);
  viewer << DGtal::CustomColors3D(DGtal::Color(0,250,0),DGtal::Color(0,250,0));
  viewer.addCylinder(p1,p2,tree.myVectSegments[1].myRadius*thickness);
  
  for (auto s : tree.myVectSegments) {
    // test if the segment is the root or its parent we do not display (already done).
    if (s.myIndex == 0 || s.myIndex == 1)
      continue;
 
    DGtal::Z3i::RealPoint distal = s.myCoordinate;
    DGtal::Z3i::RealPoint proxital = tree.myVectSegments[tree.myVectParent[s.myIndex]].myCoordinate;
    viewer << DGtal::CustomColors3D(DGtal::Color(250,0,0),DGtal::Color(250,0,0));
    viewer.addBall(distal,tree.myVectSegments[s.myIndex].myRadius);
    viewer << DGtal::CustomColors3D(DGtal::Color(0,250,0),DGtal::Color(0,250,0));
    //viewer.addBall(distal,1);
    viewer.addCylinder(distal,proxital,tree.myVectSegments[s.myIndex].myRadius*thickness);
    //std::cout<<"r="<<tree.myVectSegments[s.myIndex].myRadius<<std::endl;
    i++;
  }
  /*
  //Display Sphere domaine
  viewer << DGtal::CustomColors3D(DGtal::Color(0,0,250,10),DGtal::Color(0,0,250,10));
  viewer.addBall(DGtal::Z3i::RealPoint(0,0,0),my_rPerf);
  */
  viewer<< MyViewer::updateDisplay;
  return application.exec();

  return EXIT_SUCCESS;
}
