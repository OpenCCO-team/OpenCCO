//
//  TestGeneration.cpp
//  
//
//  Created by Dylan on 16/07/2020.
//


#include <iostream>
#include <stdio.h>
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
        if(cRand.GetLength(cRand.getNearestSegment(NewPoint)) !=0 ){
            
        cRand.addSegmentFromPoint(NewPoint);
            
        cRand.myKTerm++;
        cRand.updateScale(cRand.myRsupp);

        }
        
    }
}
        cout << " a la sortie, le segment 0 et 1 ont pour radius" << cRand.myVectParent[0] << " and " << cRand.myVectParent[1] << endl;

    

    cRand.exportDisplay("cooAlgo_testGene.svg");
    DGtal::trace.info() << "Parsing done..." << std::endl;
    return EXIT_SUCCESS;
}


