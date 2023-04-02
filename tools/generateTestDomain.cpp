/**
 *  GenerateTestDomain program (used in OpenCCO implementation) 
 *  Copyright (C) 2023 B. Kerautret;  Phuc Ngo, N. Passat H. Talbot and C. Jaquet
 *
 *  This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 **/

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

