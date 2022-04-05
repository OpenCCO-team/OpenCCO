#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <chrono>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/io/Display3D.h"
#include "DGtal/io/viewers/Viewer3D.h"

#include "CoronaryArteryTree.h"
#include "geomhelpers.h"

/**
 * @brief read seed file
 * @param filename
 * @return vector of pair of point and its corresponding radius
 */
std::vector<std::pair<CoronaryArteryTree::Point3D, double> >
readSeed(std::string filename) {
  std::vector<std::pair<CoronaryArteryTree::Point3D, double> > vecSeeds;
  double x, y, r;
  std::ifstream myfile (filename);
  if (myfile.is_open()) {
    while ( myfile >> x >> y >> r )
      vecSeeds.push_back(std::make_pair(CoronaryArteryTree::Point3D(x,y), r));
    myfile.close();
  }
  else std::cout << "Unable to open file"<<std::endl;
  assert(vecSeeds.size() != 0);
  return vecSeeds;
}

/**
 * @brief read seed file
 * @param filename
 * @return vector of pair of point and its corresponding radius
 */

CoronaryArteryTree
testAutoGen(double aPerf, int nbTerm) {
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test random adds with distance constraint");
  srand (time(NULL));
  double rRoot = 1.0;//10.0/nbTerm;
  std::string filename;
  
  CoronaryArteryTree cTree (aPerf, nbTerm, rRoot);
  
  bool isOK = false;
  unsigned int nbSeed = cTree.my_NTerm;
  for (unsigned int i = 1; i < nbSeed; i++) {
    DGtal::trace.progressBar(i, nbSeed);
    int nbSol = 0, itOpt = 0;
    CoronaryArteryTree cTreeOpt = cTree;
    double volOpt = -1.0, vol = 0.0;
    while (nbSol==0) {
      CoronaryArteryTree::Point3D pt = cTree.generateNewLocation(100);
      std::vector<unsigned int> vecN = cTree.getN_NearestSegments(pt,cTree.myNumNeighbor);
      for(size_t it=0; it<vecN.size(); it++) {
        //if(!cTree.isIntersecting(pt, cTree.findBarycenter(pt, vecN.at(it)),vecN.at(it),n))
        if(!cTree.isIntersecting(pt, cTree.findBarycenter(pt, vecN.at(it)),vecN.at(it),cTree.myNumNeighbor, 2*cTree.myVectSegments[vecN.at(it)].myRadius)) {
          CoronaryArteryTree cTree1 = cTree;
          isOK = cTree1.isAddable(pt,vecN.at(it), 100, 0.01, cTree1.myNumNeighbor);
          if(isOK) {
            vol = cTree1.computeTotalVolume(1);
            if(volOpt<0.0) {
              volOpt = vol;
              cTreeOpt = cTree1;
              itOpt = it;
            }
            else {
              if(volOpt>vol) {
                volOpt = vol;
                cTreeOpt = cTree1;
                itOpt = it;
              }
            }
            nbSol++;
          }
        }
      }
    }
    cTree = cTreeOpt;
    cTree.updateLengthFactor();
    cTree.updateResistanceFromRoot();
    cTree.updateRootRadius();
  }
  std::cout<<"====> Aperf="<<cTree.myRsupp*cTree.myRsupp*cTree.my_NTerm*M_PI<<" == "<<aPerf<<std::endl;

  //filename = "testCCO_"+std::to_string(nbTerm)+".eps";
  //cTree.exportBoardDisplay(filename.c_str(), 1.0);
  //cTree.myBoard.clear();
  return cTree;
}

CoronaryArteryTree
testFixedInit_AutoGen(double aPerf, int nbTerm, std::vector<CoronaryArteryTree::Point3D> constraintSeeds) {
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test random adds with distance constraint");
  srand (time(NULL));
  double rRoot = 1.0;//10.0/nbTerm;
  std::string filename;
  
  //CoronaryArteryTree cTree (aPerf, nbTerm, rRoot);
  CoronaryArteryTree cTree (aPerf, nbTerm, constraintSeeds); //init with an array of seeds
  
  //Generate randomly other seeds
  bool isOK = false;
  unsigned int nbSeed = cTree.my_NTerm;
  for (unsigned int i = 1; i < nbSeed; i++) {
    DGtal::trace.progressBar(i, nbSeed);
    int nbSol = 0, itOpt = 0;
    CoronaryArteryTree cTreeOpt = cTree;
    double volOpt = -1.0, vol = 0.0;
    while (nbSol==0) {
      CoronaryArteryTree::Point3D pt = cTree.generateNewLocation(100);
      std::vector<unsigned int> vecN = cTree.getN_NearestSegments(pt,cTree.myNumNeighbor);
      for(size_t it=0; it<vecN.size(); it++) {
        //if(!cTree.isIntersecting(pt, cTree.findBarycenter(pt, vecN.at(it)),vecN.at(it),n))
        if(!cTree.isIntersecting(pt, cTree.findBarycenter(pt, vecN.at(it)),vecN.at(it),cTree.myNumNeighbor, 2*cTree.myVectSegments[vecN.at(it)].myRadius)) {
          CoronaryArteryTree cTree1 = cTree;
          isOK = cTree1.isAddable(pt,vecN.at(it), 100, 0.01, cTree1.myNumNeighbor);
          if(isOK) {
            vol = cTree1.computeTotalVolume(1);
            if(volOpt<0.0) {
              volOpt = vol;
              cTreeOpt = cTree1;
              itOpt = it;
            }
            else {
              if(volOpt>vol) {
                volOpt = vol;
                cTreeOpt = cTree1;
                itOpt = it;
              }
            }
            nbSol++;
          }
        }
      }
    }
    cTree = cTreeOpt;
    cTree.updateLengthFactor();
    cTree.updateResistanceFromRoot();
    cTree.updateRootRadius();
  }
  std::cout<<"====> Aperf="<<cTree.myRsupp*cTree.myRsupp*cTree.my_NTerm*M_PI<<" == "<<aPerf<<std::endl;
  
  //filename = "testCCO_"+std::to_string(nbTerm)+".eps";
  //cTree.exportBoardDisplay(filename.c_str(), 1.0);
  //cTree.myBoard.clear();
  return cTree;
}

/**
 * @brief main function call
 *
 */
int main(int argc, char** argv)
{
  QApplication application(argc,argv);
  
  clock_t start, end;
  
  start = clock();
  //1000 => Execution time: 129.17274900 sec
  //2000 => Execution time: 478.48590200 sec
  //3000 => Execution time: 1023.94746700 sec
  //4000 => Execution time: 1896.94450700 sec
  //5000 => Execution time: 3435.08630500 sec
  double aPerf = 200000;
  double nTerm = 500;
  //CoronaryArteryTree tree = testAutoGen(aPerf, nTerm);
  
  double my_rPerf = pow(3.0*aPerf/(4.0*M_PI),1.0/3.0);//2D: sqrt(aPerf/M_PI);
  std::vector<CoronaryArteryTree::Point3D> constraintSeeds;
  constraintSeeds.push_back(CoronaryArteryTree::Point3D(0, my_rPerf, 0));//Root
  int n=4;
  for(int i=1; i<n; i++) {
    double angle = i*2*M_PI/n + M_PI/2.0;
    std::cout<<"iter"<<i<<", angle="<<angle<<std::endl;
    //constraintSeeds.push_back(CoronaryArteryTree::Point3D(cos(angle)*my_rPerf, sin(angle)*my_rPerf));
    constraintSeeds.push_back(CoronaryArteryTree::Point3D(sqrt(3.0)*cos(angle)*my_rPerf/2.0, -1.0/2.0, sqrt(3.0)*sin(angle)*my_rPerf)/2.0);
  }
  CoronaryArteryTree tree = testFixedInit_AutoGen(aPerf, nTerm,constraintSeeds);
  
  end = clock();
  printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);

  typedef DGtal::Viewer3D<> MyViewer;
  MyViewer viewer;
  viewer.show();
  unsigned int i = 0;
  double thickness = 1;
  viewer << DGtal::CustomColors3D(DGtal::Color(0,0,250),DGtal::Color(0,0,250));
  CoronaryArteryTree::Point3D p1 = tree.myVectSegments[1].myCoordinate;
  CoronaryArteryTree::Point3D p2 = tree.myVectSegments[tree.myVectParent[tree.myVectSegments[1].myIndex]].myCoordinate;
  viewer.addBall(p2,tree.myVectSegments[1].myRadius);
  viewer << DGtal::CustomColors3D(DGtal::Color(0,250,0),DGtal::Color(0,250,0));
  viewer.addCylinder(p1,p2,tree.myVectSegments[1].myRadius*thickness);
  
  for (auto s : tree.myVectSegments) {
    // test if the segment is the root or its parent we do not display (already done).
    if (s.myIndex == 0 || s.myIndex == 1)
      continue;
 
    CoronaryArteryTree::Point3D distal = s.myCoordinate;
    CoronaryArteryTree::Point3D proxital = tree.myVectSegments[tree.myVectParent[s.myIndex]].myCoordinate;
    viewer << DGtal::CustomColors3D(DGtal::Color(250,0,0),DGtal::Color(250,0,0));
    viewer.addBall(distal,tree.myVectSegments[s.myIndex].myRadius);
    viewer << DGtal::CustomColors3D(DGtal::Color(0,250,0),DGtal::Color(0,250,0));
    //viewer.addBall(distal,1);
    viewer.addCylinder(distal,proxital,tree.myVectSegments[s.myIndex].myRadius*thickness);
    std::cout<<"r="<<tree.myVectSegments[s.myIndex].myRadius<<std::endl;
    i++;
  }
  /*
  //Display Sphere domaine
  viewer << DGtal::CustomColors3D(DGtal::Color(0,0,250,10),DGtal::Color(0,0,250,10));
  viewer.addBall(CoronaryArteryTree::Point3D(0,0,0),my_rPerf);
  */
  viewer<< MyViewer::updateDisplay;
  return application.exec();
  
  return EXIT_SUCCESS;
}
