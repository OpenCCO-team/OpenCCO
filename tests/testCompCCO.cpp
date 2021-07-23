#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <chrono>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "CoronaryArteryTree.h"
#include "GeomHelpers.h"

/**
 * @brief reads seed file from results obtained by Clara and Hugues's code
 * @param filename : input filename
 * @return vector of pair of point and its corresponding radius
 */
std::vector<std::pair<DGtal::Z2i::RealPoint, double> >
readSeed(std::string filename) {
  std::vector<std::pair<DGtal::Z2i::RealPoint, double> > vecSeeds;
  double x, y, r;
  std::ifstream myfile (filename);
  if (myfile.is_open()) {
       while ( myfile >> x >> y >> r )
      vecSeeds.push_back(std::make_pair(DGtal::Z2i::RealPoint(x,y), r));
    myfile.close();
   
  }
  else std::cout << "Unable to open file"<<std::endl;
  assert(vecSeeds.size() != 0);
  return vecSeeds;
}

/**
 * @brief compares results obtained by Clara and Hugues's code and ours implementation
 * @param NTerm number of terminal
 * @param seed seed used for the generation of terminals
 * @return vector of pair of point and its corresponding radius
 */

void
testCompareResult(int NTerm, int seed)
{
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test adds fixed terminal points from file");
  
  std::string dir = "../Data/Nt" + std::to_string(NTerm) + "_kt10_s" + std::to_string(seed) + "_M301/";
  std::string prefix = "NoCom_Nt" + std::to_string(NTerm) + "_s" + std::to_string(seed) + "_M301_";
  std::string fileDistal = dir + prefix + "distal.txt";
  std::string fileProximal = dir + prefix + "proximal.txt";
  std::string fileSeeds = dir + prefix + "TerminalSeeds.txt";
  std::vector<std::pair<DGtal::Z2i::RealPoint, double> > vecSeed = readSeed(fileSeeds);
  
  double radius = 50;
  double aPerf = M_PI*radius*radius;//2000000;
  int nbTerm = vecSeed.size();//10;
  double rRoot = 1.0;//1.0/nbTerm;
  double r = sqrt(aPerf/(nbTerm*M_PI));
  DGtal::Z2i::RealPoint pRoot(2*radius,2*radius);
  DGtal::Z2i::RealPoint pCenter(2*radius,radius);
  //Test constructors
  //CoronaryArteryTree cTree (aPerf, nbTerm, rRoot);
  //CoronaryArteryTree cTree (pRoot, aPerf, nbTerm, rRoot);
  std::string filename;
  DGtal::Z2i::RealPoint pTerm = vecSeed[0].first;
  CoronaryArteryTree cTree (pCenter, pRoot, pTerm, nbTerm, rRoot);
  
  std::cout<<"Vol : "<<cTree.computeTotalVolume(1)<<std::endl;
  
  bool isOK = false;
  unsigned int nbSeed = cTree.my_NTerm;
  for (unsigned int i = 1; i < nbSeed; i++) {
    DGtal::trace.progressBar(i, nbSeed);
    int nbSol = 0, itOpt = 0;
    CoronaryArteryTree cTreeOpt = cTree;
    double volOpt = -1.0, vol = cTree.computeTotalVolume();
    std::cout<<"Vol in ("<<i<<"): "<< vol <<std::endl;
    while (nbSol==0) {
      CoronaryArteryTree::Point2D pt = vecSeed[i].first; //cTree.generateNewLocation(100);
      std::vector<unsigned int> vecN = cTree.getN_NearestSegments(pt,cTree.myNumNeighbor);
      for(size_t it=0; it<vecN.size(); it++) {
        if(!cTree.isIntersecting(pt, cTree.findBarycenter(pt, vecN.at(it)),vecN.at(it),cTree.myNumNeighbor, cTree.myVectSegments[vecN.at(it)].myRadius)) {
          CoronaryArteryTree cTree1 = cTree;
          isOK = cTree1.isAddable(pt,vecN.at(it), 100, 0.01, cTree1.myNumNeighbor);
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
    cTree.updateLengthFactor();
    cTree.updateResistanceFromRoot();
    cTree.updateRootRadius();
    vol = cTree.computeTotalVolume();
    std::cout<<"Vol out ("<<i<<"): "<< vol <<std::endl;
  }
  //std::cout<<"====> Aperf="<<cTree.myRsupp*cTree.myRsupp*cTree.my_NTerm*M_PI<<" == "<<aPerf<<std::endl;

  //Draw CCO result
  std::vector<std::pair<DGtal::Z2i::RealPoint, double> > vecCCO_res1 = readSeed(fileDistal);
  std::vector<std::pair<DGtal::Z2i::RealPoint, double> > vecCCO_res2 = readSeed(fileProximal);
  //from Clara and Hugues's code
  for(size_t it=0; it<vecCCO_res1.size(); it++) {
      DGtal::Z2i::RealPoint p1 = vecCCO_res1.at(it).first;
      DGtal::Z2i::RealPoint p2 = vecCCO_res2.at(it).first;
      double r = vecCCO_res1.at(it).second;
      cTree.myBoard.setPenColor(DGtal::Color::Black);
    cTree.myBoard.fillCircle(p2[0], p2[1], 20*r/57.5, 1);
    cTree.myBoard.setPenColor(DGtal::Color::Green);
    cTree.myBoard.setLineWidth(20*r);
    cTree.myBoard.drawLine(p1[0], p1[1], p2[0], p2[1], 2);
  }
  filename = "testCCO_Nt" + std::to_string(NTerm) + "_s" + std::to_string(seed) +".eps";
  //our result
  cTree.exportBoardDisplay(filename.c_str(), 1.0, true, false);
  cTree.myBoard.clear();

  
}
/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{
  clock_t start, end;
  int Nt = 10; //10 20 30 40 50 60
  int seed = 42;//42 420 25 250 90 201 15 215
  start = clock();
  testCompareResult(Nt, seed);
  end = clock();
  printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
  
  return EXIT_SUCCESS;
}
