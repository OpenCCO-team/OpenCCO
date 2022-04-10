#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

#include "CoronaryArteryTree.h"
#include "GeomHelpers.h"
#include "ConstructionHelpers.h"

/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{
  std::string resource_dir = SAMPLE_DIR;
  
  ///-----------------------------------------------------------------------------------------------
  /// Test projection distances
  DGtal::trace.beginBlock("Testing distance map from 3D image domain");
  std::cout << "ressource:"  << resource_dir << std::endl;
  std::stringstream ss;
  ss << resource_dir <<"shape3D.vol";
  CoronaryArteryTree<2>::Image img = DGtal::GenericReader<CoronaryArteryTree<2>::Image>::import(ss.str());

  DGtal::trace.beginBlock("Testing computation of ditance map..");
  typedef typename DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> ImageDouble;
  ImageDouble imgD = GeomHelpers::getImageDistance2D<CoronaryArteryTree<2>::Image, ImageDouble >(img);
  DGtal::Z2i::Point pInt(313, 201);
  DGtal::Z2i::Point pExt(20, 20);
  
  
  
  DGtal::trace.info() << " Distance map of point int "  << pInt << "distance:" <<
  (int) imgD(pInt) << " sould be > 10" << (imgD(pInt) > 10 ? " OK": "KO") << std::endl;
  
  DGtal::trace.info() << " Distance map of point int "  << pExt << "distance:" <<
  (int) imgD(pExt) << " sould be 0" << ( imgD(pExt) == 0 ? " OK": "KO") << std::endl;

  DGtal::trace.info() << "Export distance map...";
  DGtal::trace.info() << "[Done]" << std::endl;
  imgD >> "distanceMap.pgm" ;
  DGtal::trace.endBlock();

  

  //  CoronaryArteryTree<2>::Image img = DGtal::GenericReader<CoronaryArteryTree<2>::Image>::import(ss.str());
  //  bool checkDomInter = GeomHelpers::checkNoIntersectDomain(img, 128, DGtal::Z2i::Point(264,196), DGtal::Z2i::Point(438,225));
  //  DGtal::trace.info() << "Test intersection: 264 196 and 438 225 "
  //  << " distance (should be true) :" << checkDomInter <<  ( checkDomInter ? " OK": " KO")  << std::endl;
  //  bool checkDomInter2 = GeomHelpers::checkNoIntersectDomain(img, 128, DGtal::Z2i::Point(264,196), DGtal::Z2i::Point(2,56));
  //  DGtal::trace.info() << "Test intersection: 264 196 and 2 56 "
  //  << " distance (should be false) :" << checkDomInter2 <<  ( !checkDomInter2 ? " OK": " KO")  << std::endl;
  //  DGtal::trace.endBlock();
  //

//  DGtal::trace.beginBlock("Testing generate circle points.");
//  DGtal::Z2i::DigitalSet vPt =  GeomHelpers::pointsOnCircle(DGtal::Z2i::Point(3,3), 4);
//  DGtal::trace.info() << "Nb points :" << vPt.size() <<
//  " Should be more than 0 "<< ( vPt.size() > 0 ? " OK": "KO") <<   std::endl;
//  DGtal::trace.info() << "Should contain P(7,3) and P(-1,3) :"
//  << (vPt.find(DGtal::Z2i::Point(7,3)) != vPt.end()  && vPt.find(DGtal::Z2i::Point(-1,3)) != vPt.end() ? " OK": "KO") << std::endl;
//
//
//  DGtal::trace.endBlock();

  
  
  
  return EXIT_SUCCESS;
}






