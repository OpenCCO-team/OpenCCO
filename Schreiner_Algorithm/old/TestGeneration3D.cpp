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
#include "CoronaryArteryTree3D.h"

#include "DGtal/base/Common.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/boards/Board3D.h"

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
//
{
//
//
    CLI::App app;
    unsigned int nbSegments {0};
    std::string outputFileName = "result.svg";
//
//
    // BEGIN parse command line using CLI ----------------------------------------------
    app.description("Implementation of the CCO algorithm (Constrained Construction Optimisation)");
    app.add_option("-n,--nbSegments", nbSegments, "Nb segments to be generated")
      -> required();
    app.add_option("-o,--outputFileName", outputFileName, "output file name (svg).", true);

//
//

    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------------


srand (time(NULL));
CoronaryArteryTree3D cRand(DGtal::Z3i::RealPoint(0, 25, 0), 2000, 1);
cRand.addFirstSegment(DGtal::Z3i::RealPoint(0, 0, 0));
cRand.Nterm=nbSegments;
//
// The Coronary artery tree is created
DGtal::trace.beginBlock("Parameters printing");
DGtal::trace.info() << cRand;
while(cRand.myKTerm < nbSegments)
{
//
    double rho=( rand()/double(RAND_MAX))*cRand.myRsupp;
    double theta =(rand()/double(RAND_MAX))*PI;
    double phi =(rand()/double(RAND_MAX))*2*PI;

    double x = rho*sin(theta)*cos(phi);
    double y = rho*sin(theta)*sin(phi);
    double z = rho*cos(theta);
    
    DGtal::Z3i::RealPoint NewPoint(x,y,z);
    if (isInsideSphere(DGtal::Z3i::RealPoint(0, 0, 0), DGtal::Z3i::RealPoint(x, y, z),  cRand.myRsupp)){
        cout << cRand.GetLength(cRand.getNearestSegment(NewPoint)) << endl;
        if(cRand.GetLength(cRand.getNearestSegment(NewPoint)) !=0 ){

        cRand.addSegmentFromPoint(NewPoint);

        cRand.myKTerm++;
        cRand.updateScale(cRand.myRsupp);

        }

    }
//    cout << "kterm = " << cRand.myKTerm<< endl;
}
//        cout << " a la sortie, le segment 0 et 1 ont pour radius" << cRand.myVectParent[0] << " and " << cRand.myVectParent[1] << endl;
//
//
//
//    cRand.exportDisplay("cooAlgo_testGene3D.svg");
    cRand.CreateVertexFile("TestVertex.txt");
    cRand.CreateRadiusFile("TestRadius.txt");
    cRand.CreateEdgesFile("TestEdges.txt");

    DGtal::trace.info() << "Parsing done..." << std::endl;
    return EXIT_SUCCESS;
}
//
//


// generique reader a prendre. Pas d'exemple propre.
