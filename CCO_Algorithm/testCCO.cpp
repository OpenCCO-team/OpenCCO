#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <chrono>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "CoronaryArteryTree.h"
#include "geomhelpers.h"

/**
 * @brief read seed file
 * @param filename
 * @return vector of pair of point and its corresponding radius
 */
std::vector<std::pair<DGtal::Z2i::RealPoint, double> > readSeed(std::string filename) {
  std::vector<std::pair<DGtal::Z2i::RealPoint, double> > vecSeeds;
  double x, y, r;
  std::ifstream myfile (filename);
  if (myfile.is_open()) {
    while ( myfile >> x >> y >> r )
      vecSeeds.push_back(std::make_pair(DGtal::Z2i::RealPoint(x,y), r));
    myfile.close();
  }
  else std::cout << "Unable to open file"<<std::endl;

  return vecSeeds;
}

void testAutoGen(double aPerf, int nbTerm) {
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test random adds with distance constraint");
  srand (time(NULL));
  double rRoot = 10.0/nbTerm;
  
  CoronaryArteryTree cTree (aPerf, nbTerm, rRoot);
  
  unsigned int n = 20;
  bool isOK = false;
  std::string filename;
  unsigned int nbSeed = cTree.my_NTerm - 1;
  for (unsigned int i = 0; i < nbSeed; i++) {
    DGtal::trace.progressBar(i, nbSeed);
    int nbSol = 0, itOpt = 0;
    CoronaryArteryTree cTreeOpt = cTree;
    double volOpt = -1.0, vol = 0.0;
    while (nbSol==0) {
      CoronaryArteryTree::Point2D pt = cTree.generateNewLocation(100);
      std::vector<unsigned int> vecN = cTree.getN_NearestSegments(pt,n);
      for(size_t it=0; it<vecN.size(); it++) {
        if(!cTree.isIntersecting(pt, cTree.FindBarycenter(pt, vecN.at(it)),vecN.at(it),n)) {
          CoronaryArteryTree cTree1 = cTree;
          isOK = cTree1.isAddable(pt,vecN.at(it), 10, n);
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
  }
  std::cout<<"====> Aperf="<<cTree.myRsupp*cTree.myRsupp*cTree.my_NTerm*M_PI<<" == "<<aPerf<<std::endl;

  filename = "testCCO_"+std::to_string(nbTerm)+".eps";
  cTree.exportBoardDisplay(filename.c_str(), 5.0, true);
  cTree.myBoard.clear();
}

void testAutoGen2(double aPerf, int nbTerm) {
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test random adds with distance constraint");
  srand (time(NULL));
  double rRoot = 10.0/nbTerm;
  std::string filename;
  
  CoronaryArteryTree cTree (aPerf, nbTerm, rRoot);
  
  unsigned int n = 20;
  bool isOK = false;
  unsigned int nbSeed = cTree.my_NTerm;
  for (unsigned int i = 1; i < nbSeed; i++) {
    DGtal::trace.progressBar(i, nbSeed);
    int nbSol = 0, itOpt = 0;
    CoronaryArteryTree cTreeOpt = cTree;
    double volOpt = -1.0, vol = 0.0;
    while (nbSol==0) {
      CoronaryArteryTree::Point2D pt = cTree.generateNewLocation(100);
      std::vector<unsigned int> vecN = cTree.getN_NearestSegments(pt,n);
      for(size_t it=0; it<vecN.size(); it++) {
        if(!cTree.isIntersecting(pt, cTree.FindBarycenter(pt, vecN.at(it)),vecN.at(it),n))
        {
          CoronaryArteryTree cTree1 = cTree;
          isOK = cTree1.isAddable(pt,vecN.at(it), 100, 0.01, n);
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
  }
  std::cout<<"====> Aperf="<<cTree.myRsupp*cTree.myRsupp*cTree.my_NTerm*M_PI<<" == "<<aPerf<<std::endl;

  filename = "testCCO_"+std::to_string(nbTerm)+".eps";
  cTree.exportBoardDisplay(filename.c_str(), 5.0, true);
  cTree.myBoard.clear();
}

void testFixedSeeds(double radius, std::vector<std::pair<DGtal::Z2i::RealPoint, double> > vecSeed) {
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test adds fixed terminal points from file");
  double rRoot = 1.0;
  double aPerf = M_PI*radius*radius;//2000000;
  int nbTerm = vecSeed.size();//10;
  DGtal::Z2i::RealPoint pRoot(2*radius,2*radius);
  DGtal::Z2i::RealPoint pCenter(2*radius,radius);
  //Test constructors
  std::string filename;
  DGtal::Z2i::RealPoint pTerm = vecSeed[0].first;
  CoronaryArteryTree cTree (pCenter, pRoot, pTerm, aPerf, nbTerm, rRoot);
  /*
  for(size_t it=1; it<vecSeed.size(); it++) {
    CoronaryArteryTree::Point2D pt = vecSeed[it].first;
    //cTree.addSegmentFromPointWithBarycenter(pt);
    cTree.addSegmentFromPoint(pt);
  }
  */
  std::cout<<"Tree vol: "<<cTree.computeTotalVolume(1)<<std::endl;
  
  unsigned int n = 20;
  bool isOK = false;
  unsigned int nbSeed = cTree.my_NTerm;
  for (unsigned int i = 1; i < nbSeed; i++) {
    DGtal::trace.progressBar(i, nbSeed);
    int nbSol = 0, itOpt = 0;
    CoronaryArteryTree cTreeOpt = cTree;
    double volOpt = -1.0, vol = 0.0;
    while (nbSol==0) {
      CoronaryArteryTree::Point2D pt = vecSeed[i].first; //cTree.generateNewLocation(100);
      std::vector<unsigned int> vecN = cTree.getN_NearestSegments(pt,n);
      for(size_t it=0; it<vecN.size(); it++) {
        if(!cTree.isIntersecting(pt, cTree.FindBarycenter(pt, vecN.at(it)),vecN.at(it),n))
        {
          CoronaryArteryTree cTree1 = cTree;
          isOK = cTree1.isAddable(pt,vecN.at(it), 100, n);
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
  }
  std::cout<<"====> Aperf="<<cTree.myRsupp*cTree.myRsupp*cTree.my_NTerm*M_PI<<" == "<<aPerf<<std::endl;
  
  //Draw CCO result
  std::vector<std::pair<DGtal::Z2i::RealPoint, double> > vecCCO_res1 = readSeed("../Data/NoCom_Nt10_s420_M301_distal.txt");//NoCom_Nt10_s420_M301_distal InterTree_Nt10_kt2_s420_M301_distal
  std::vector<std::pair<DGtal::Z2i::RealPoint, double> > vecCCO_res2 = readSeed("../Data/NoCom_Nt10_s420_M301_proximal.txt");//NoCom_Nt10_s420_M301_proximal InterTree_Nt10_kt2_s420_M301_proximal
  
  for(size_t it=0; it<vecCCO_res1.size(); it++) {
    DGtal::Z2i::RealPoint p1 = vecCCO_res1.at(it).first;
    DGtal::Z2i::RealPoint p2 = vecCCO_res2.at(it).first;
    double r = vecCCO_res1.at(it).second;
    cTree.myBoard.setPenColor(DGtal::Color::Black);
    cTree.myBoard.fillCircle(p2[0], p2[1], 20*r/57.5, 1);
    cTree.myBoard.setPenColor(DGtal::Color::Green);
    cTree.myBoard.setLineWidth(20*r);
    cTree.myBoard.drawLine(p1[0], p1[1], p2[0], p2[1], 2);
    //cTree.myBoard.setPenColor(DGtal::Color::Red);
    //cTree.myBoard.drawLine(pRoot[0], pRoot[1], p1[0], p1[1], 2);
  }
  filename = "testCCO_"+std::to_string(nbTerm)+".eps";
  cTree.exportBoardDisplay(filename.c_str(), 1.0, true, false);
  cTree.myBoard.clear();
}

/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{
  std::vector<std::pair<DGtal::Z2i::RealPoint, double> > vecSeed = readSeed("../Data/Nt10_kt10_s42_M301/NoCom_Nt10_s42_M301_TerminalSeeds.txt");
  assert(vecSeed.size() != 0);
  
  clock_t start, end;
  start = clock();
  //testAutoGen2(20000, 100);
  testAutoGen2(200000, 1000);
  //testFixedSeeds(50, vecSeed);
  end = clock();
  printf ("Execution time: %0.8f sec\n", ((float) end - start)/CLOCKS_PER_SEC);
  return 0;
  
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test adds fixed terminal points from file");
  srand (time(NULL));
  
  double radius = 50;
  double aPerf = M_PI*radius*radius;//2000000;
  int nbTerm = vecSeed.size();//10;
  double rRoot = 1.0;//1.0/nbTerm;
  double r = sqrt(aPerf/(nbTerm*M_PI));
  DGtal::Z2i::RealPoint pRoot(2*radius,2*radius);
  DGtal::Z2i::RealPoint pCenter(2*radius,radius);
  //DGtal::Z2i::RealPoint pTerm(0,0);
  //Test constructors
  //CoronaryArteryTree cTree (aPerf, nbTerm, rRoot);
  //CoronaryArteryTree cTree (pRoot, aPerf, nbTerm, rRoot);
  std::string filename;
  DGtal::Z2i::RealPoint pTerm = vecSeed[0].first;
  CoronaryArteryTree cTree (pCenter, pRoot, pTerm, aPerf, nbTerm, rRoot);
  
  std::cout<<"Vol : "<<cTree.computeTotalVolume(1)<<std::endl;
  
  unsigned int n = 20;
  bool isOK = false;
  unsigned int nbSeed = cTree.my_NTerm;
  for (unsigned int i = 1; i < nbSeed; i++) {
    DGtal::trace.progressBar(i, nbSeed);
    int nbSol = 0, itOpt = 0;
    CoronaryArteryTree cTreeOpt = cTree;
    double volOpt = -1.0, vol = 0.0;
    std::cout<<"Vol in ("<<i<<"): "<<cTree.computeTotalVolume()<<std::endl;
    while (nbSol==0) {
      CoronaryArteryTree::Point2D pt = vecSeed[i].first; //cTree.generateNewLocation(100);
      std::vector<unsigned int> vecN = cTree.getN_NearestSegments(pt,n);
      for(size_t it=0; it<vecN.size(); it++) {
        if(!cTree.isIntersecting(pt, cTree.FindBarycenter(pt, vecN.at(it)),vecN.at(it),n, cTree.myVectSegments[vecN.at(it)].myRadius)) {
          CoronaryArteryTree cTree1 = cTree;
          isOK = cTree1.isAddable(pt,vecN.at(it), 100, 0.01, n);
          if(isOK) {
            vol = cTree1.computeTotalVolume();
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
    //std::cout<<"Vol out 1: "<<cTree.computeTotalVolume()<<std::endl;
    cTree.updateLengthFactor();
    //cTree.updateResistanceTerminal(cTree.myVectSegments.size()-1); //seg right (new segment)
    //cTree.updateResistance(cTree.myVectSegments.size()-2); //seg left
    //cTree.updateResistance(cTree.myVectSegments.size()-1, 0);
    cTree.updateResistanceFromRoot();
    //cTree.updateBeta(cTree.myVectSegments.size()-1);
    //cTree.updateBeta();
    cTree.updateRootRadius();
    //std::cout<<"Vol out 2: "<<cTree.computeTotalVolume()<<std::endl;
    std::cout<<"Vol out ("<<i<<"): "<<cTree.computeTotalVolume()<<std::endl;
  }
  std::cout<<"====> Aperf="<<cTree.myRsupp*cTree.myRsupp*cTree.my_NTerm*M_PI<<" == "<<aPerf<<std::endl;

  //Draw CCO result
  std::vector<std::pair<DGtal::Z2i::RealPoint, double> > vecCCO_res1 = readSeed("../Data/Nt10_kt10_s42_M301/InterTree_Nt10_kt10_s42_M301_distal.txt");//NoCom_Nt10_s42_M301_distal InterTree_Nt10_kt3_s42_M301_distal
  std::vector<std::pair<DGtal::Z2i::RealPoint, double> > vecCCO_res2 = readSeed("../Data/Nt10_kt10_s42_M301/InterTree_Nt10_kt10_s42_M301_proximal.txt");//NoCom_Nt10_s42_M301_proximal InterTree_Nt10_kt3_s42_M301_proximal
  
  for(size_t it=0; it<vecCCO_res1.size(); it++) {
    DGtal::Z2i::RealPoint p1 = vecCCO_res1.at(it).first;
    DGtal::Z2i::RealPoint p2 = vecCCO_res2.at(it).first;
    double r = vecCCO_res1.at(it).second;
    cTree.myBoard.setPenColor(DGtal::Color::Black);
    cTree.myBoard.fillCircle(p2[0], p2[1], 20*r/57.5, 1);
    cTree.myBoard.setPenColor(DGtal::Color::Green);
    cTree.myBoard.setLineWidth(20*r);
    cTree.myBoard.drawLine(p1[0], p1[1], p2[0], p2[1], 2);
    //cTree.myBoard.setPenColor(DGtal::Color::Red);
    //cTree.myBoard.drawLine(pRoot[0], pRoot[1], p1[0], p1[1], 2);
  }
  filename = "testCCO_"+std::to_string(nbTerm)+".eps";
  cTree.exportBoardDisplay(filename.c_str(), 1.0, true, false);
  cTree.myBoard.clear();

  return EXIT_SUCCESS;
}
