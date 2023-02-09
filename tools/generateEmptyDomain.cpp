
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "DGtal/base/Common.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/shapes/implicit/ImplicitBall.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/shapes/EuclideanShapesDecorator.h"

#include "CLI11.hpp"


using namespace std;


int
main(int argc, char** argv)

{
    typedef DGtal::ImageContainerBySTLVector < DGtal::Z3i::Domain, unsigned char> Image3D;
    DGtal::Z3i::RealPoint ptL(0,0,0);
    DGtal::Z3i::RealPoint ptU(50,50,50);
    Image3D result(DGtal::Z3i::Domain(ptL, ptU));
    for (auto p: result.domain())
    {
      if (p[0] > 5 && p[1] > 5 && p[2] > 5 && p[0] < 45 && p[1] < 45&& p[2] < 45   )
      result.setValue(p, 200);
    }
    result >> "empty.vol";
    
    return EXIT_SUCCESS;
}

