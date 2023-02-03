
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
    typedef DGtal::ImplicitBall< DGtal::Z3i::Space > ImpBall;
    typedef DGtal::ImageContainerBySTLVector < DGtal::Z3i::Domain, unsigned char> Image3D;

    typedef DGtal::EuclideanShapesCSG< ImpBall, ImpBall > Addition;
    typedef DGtal::GaussDigitizer< DGtal::Z3i::Space, Addition> GaussDigi;

    DGtal::Z3i::RealPoint ptL(0,0,0);
    DGtal::Z3i::RealPoint ptU(200,200,200);
    DGtal::Z3i::RealPoint ptc(100,100,100);
    DGtal::Z3i::RealPoint ptc2(120,150,130);

    DGtal::Z3i::Domain dom(ptL, ptU);
    ImpBall aBall1(ptc, 50);
    ImpBall aBall2(ptc2, 50);
    Addition addShape(aBall1);
    addShape.plus(aBall2);
    
    GaussDigi digShape;
    digShape.attach( addShape );
    digShape.init(addShape.getLowerBound(), addShape.getUpperBound(),1);
   
    DGtal::Z3i::Domain domainShape = digShape.getDomain();
    DGtal::Z3i::DigitalSet aSet( domainShape );
    DGtal::Shapes<DGtal::Z3i::Domain>::digitalShaper( aSet, digShape );
    
    Image3D result(domainShape);
    
    for (auto s : aSet){
        result.setValue(s, 255);
    }
    result >> "result.vol";
    
    return EXIT_SUCCESS;
}

