//
//  TestGeneration.cpp
//
//
//  Created by Dylan on 16/07/2020.
//


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

#include "geomhelpers.h"

#include "CLI11.hpp"
#define PI 3.14159265
using namespace std;


int
main(int argc, char** argv)

{
    

    CLI::App app;
    unsigned int nbSegments {0};
    std::string outputFileName = "result.svg";


    // BEGIN parse command line using CLI ----------------------------------------------
    app.description("Implementation of the CCO algorithm (Constrained Construction Optimisation)");
    app.add_option("-n,--nbSegments", nbSegments, "Nb segments to be generated")
      -> required();
    app.add_option("-o,--outputFileName", outputFileName, "output file name (svg).", true);




    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------------
    
    
srand (time(NULL));
CoronaryArteryTree cRand (DGtal::Z2i::RealPoint(0, 25), 2000, 1);
cRand.addFirstSegment(DGtal::Z2i::RealPoint(0, 0));
cRand.Nterm=nbSegments;
    int Ntoss=0;
// The Coronary artery tree is created
DGtal::trace.beginBlock("Parameters printing");
DGtal::trace.info() << cRand;
    
    
while(cRand.myKTerm < nbSegments)
{
    
    double r=( rand()/double(RAND_MAX))*cRand.myRsupp;
    double theta =(rand()/double(RAND_MAX))*2*PI;
        
    double x= r*cos(theta);
    double y= r*sin(theta);
    DGtal::Z2i::RealPoint NewPoint(x,y);
    if (isInsideCircle(DGtal::Z2i::RealPoint(0, 0), DGtal::Z2i::RealPoint(x, y),  cRand.myRsupp)){
            DGtal::Z2i::RealPoint FinalPosition;
        int Index_Segment = cRand.FindOptimalSegment(NewPoint);
        DGtal::Z2i::RealPoint Barycenter = cRand.FindBarycenter(NewPoint, Index_Segment);

//        bool Intersection = hasIntersection( cRand.myVectSegments[Index_Segment].myCoordinate, cRand.myVectSegments[cRand.myVectParent[Index_Segment]].myCoordinate,NewPoint, Barycenter);


            
        if(cRand.compDistCriteria(DGtal::Z2i::RealPoint(x, y), Index_Segment)>cRand.myDThresold){
            cout << "kterm = " << cRand.myKTerm<< endl;
            cRand.addSegmentFromPoint(NewPoint,Index_Segment);
            
            cRand.myKTerm++;
            cRand.updateScale(cRand.myRsupp);

        }
        else
            Ntoss++;
        if(Ntoss>10){
            Ntoss =0;
            cRand.myDThresold*=0.9;
//            cout << "treshold update"<<endl;
        }
        
    }
}
//    cRand.updateRadius();
//    for (auto s : cRand.myVectSegments)
//    {
//        cout << " Generation : " << s.Generation << endl;
//        cout << " Radius = " << s.radius << endl;
//        cout << " Length = " << cRand.GetLength(s.index) << endl;
//    }
//
    cRand.exportDisplay("cooAlgo_testGeneTreshold.svg");
    DGtal::trace.info() << "Parsing done..." << std::endl;
    return EXIT_SUCCESS;
}


