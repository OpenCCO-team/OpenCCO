//
//  TestImage.cpp
//  
//
//  Created by Dylan on 09/09/2020.
//

#include <stdio.h>



#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "CoronaryArteryTree.h"

#include "DGtal/base/Common.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/helpers/StdDefs.h"
#include "Board/Shapes.h"
#include "DGtal/images/ArrayImageIterator.h"
#include "DGtal/shapes/Shapes.h"
#include <vector>
#include "DGtal/io/readers/PPMReader.h"
#include "DGtal/io/readers/PGMReader.h"
#include "geomhelpers.h"

#include "CLI11.hpp"
#define PI 3.14159265
using namespace std;
using namespace DGtal::Z2i;



int
main(int argc, char** argv)

{
  
  
  CLI::App app;
  unsigned int nbSegments {0};
  std::string outputFileName = "result.svg";
  std::string inputFileName = "Shape.pgm";
  
  // BEGIN parse command line using CLI ----------------------------------------------
  app.description("Implementation of the CCO algorithm (Constrained Construction Optimisation)");
  app.add_option("-n,--nbSegments", nbSegments, "Nb segments to be generated")
  -> required();
  app.add_option("-o,--outputFileName", outputFileName, "output file name (svg).", true);
  app.add_option("-i,--inputFileName", inputFileName, "input file name (pgm).", true);
  
  
  
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------
  
  
  srand (time(NULL));
  CoronaryArteryTree cRand (DGtal::Z2i::RealPoint(0, 25), 200000, 1);
  cRand.addFirstSegment(DGtal::Z2i::RealPoint(0, 0));
  cRand.Nterm=nbSegments;
  int Ntoss=0;
  // The Coronary artery tree is created
  DGtal::trace.beginBlock("Parameters printing");
  DGtal::trace.info() << cRand;
  
  typedef DGtal::ImageSelector < DGtal::Z2i::Domain, unsigned int>::Type Image;
  Image image = DGtal::PGMReader<Image>::importPGM( inputFileName );
  double xmax,ymax;
  xmax = cRand.FindXmax(image.domain().upperBound()[0],image.domain().upperBound()[1]);
  ymax = cRand.FindYmax(image.domain().upperBound()[0],image.domain().upperBound()[1]);
  
  
  while(cRand.myKTerm < nbSegments)
  {
    
    double r =( rand()/double(RAND_MAX))*cRand.myRsupp;
    double theta =(rand()/double(RAND_MAX))*2*PI;
    
    double x = r*cos(theta);
    double y = r*sin(theta);
    DGtal::Z2i::RealPoint NewPoint(x,y);
    
    if(x <= xmax && y <= ymax && x >= -xmax && y >= -ymax ){
      int ximage;
      int yimage;
      ximage = -(x-xmax)*(image.domain().upperBound()[0]/2)/xmax; // This has a value between 0 and upperBound
      yimage = -(y-ymax)*(image.domain().upperBound()[1]/2)/ymax;
      
      if(image({ximage,yimage}) == 0 )
      {
        int Index_Segment = cRand.FindOptimalSegment(NewPoint);
        
        if (Index_Segment == -1) {
          Ntoss++;
          
          if(Ntoss>10){
            Ntoss =0;
            cRand.myDThresold*=0.9;
          }
          continue;
        }
        DGtal::Z2i::RealPoint Barycenter = cRand.FindBarycenter(NewPoint, Index_Segment);
        //               bool Intersection = hasIntersection( cRand.myVectSegments[Index_Segment].myCoordinate, cRand.myVectSegments[cRand.myVectParent[Index_Segment]].myCoordinate,NewPoint, Barycenter);
        
        bool Intersection = cRand.hasIntersections( cRand.myVectSegments[Index_Segment], NewPoint);
        
        
        if( Index_Segment!=-1  && cRand.compDistCriteria(DGtal::Z2i::RealPoint(x, y), Index_Segment) > cRand.myDThresold && !Intersection ){
          cout << "kterm = " << cRand.myKTerm<< endl;
          cRand.addSegmentFromPoint(NewPoint,Index_Segment);
          
          cRand.myKTerm++;
          cRand.updateScale(cRand.myRsupp);
          xmax = cRand.FindXmax(image.domain().upperBound()[0],image.domain().upperBound()[1]);
          ymax = cRand.FindYmax(image.domain().upperBound()[0],image.domain().upperBound()[1]);
        }
        else
          Ntoss++;
        
        if(Ntoss>10){
          Ntoss =0;
          cRand.myDThresold*=0.9;
        }
        
        
      }
      
    }
  }
  
  DGtal::trace.info() << "Image 3D = "<<image<<std::endl;
  //
  //    for (auto s : cRand.myVectSegments)
  //    {
  //        cout << " Generation : " << s.Generation << endl;
  //        cout << " Radius = " << s.radius << endl;
  //        cout << " Length = " << cRand.GetLength(s.index) << endl;
  //    }
  cRand.exportDisplay("ccoAlgo_testGeneImage.svg");
  DGtal::trace.info() << "Parsing done..." << std::endl;
  return EXIT_SUCCESS;
  
  
}
