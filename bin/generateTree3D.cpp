#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

#include <iostream>
#include <fstream>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/Mesh.h"
#include "DGtal/io/writers/MeshWriter.h"

#include "CLI11.hpp"

#include "CoronaryArteryTree.h"
#include "GeomHelpers.h"
#include "ExpandTreeHelpers.h"
#include "XmlHelpers.h"

#include "DomainController.h"

#ifdef WITH_VISU3D_QGLVIEWER
#include "DGtal/io/viewers/Viewer3D.h"
#endif
/**
 * @brief main function call
 *
 */

/**
 * Function to construct the tree with the help ConstructionHelpers by using a domain of reconstruction defined fram an image (ImageMaskDomainCtrl)l
 */
template<typename TTree>
void
constructTreeMaskDomain(TTree &aTree, int distSearchRoot,
                       bool verbose)
{
    clock_t start, end;
    start = clock();
    ExpandTreeHelpers::initFirtElemTree(aTree, distSearchRoot);
    ExpandTreeHelpers::expandTree(aTree, verbose);
    end = clock();
    printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
}


template<typename TTree>
void
constructTreeImplicitDomain(TTree &aTree, bool verbose)
{
    clock_t start, end;
    start = clock();
    ExpandTreeHelpers::expandTree(aTree, verbose);
    end = clock();
    printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
}




template <typename TTree>
void exportResultingMesh(const TTree &tree, std::string outName)
{
    unsigned int i = 0;
    double thickness = 1;
    // Export 3D mesh of the tree
    DGtal::Mesh<DGtal::Z3i::RealPoint> aMesh(true);
    for (auto s : tree.myVectSegments) {
      // test if the segment is the root or its parent we do not display (already done).
      if (s.myIndex == 0 || s.myIndex == 1)
        continue;
      DGtal::Z3i::RealPoint distal = s.myCoordinate;
      DGtal::Z3i::RealPoint proxital = tree.myVectSegments[tree.myVectParent[s.myIndex]].myCoordinate;
      auto v = {distal, proxital};
      DGtal::Mesh<DGtal::Z3i::RealPoint>::createTubularMesh(aMesh, v, tree.myVectSegments[s.myIndex].myRadius*thickness, 0.05);
      i++;
    }
    aMesh >> outName;
}

template <typename TTree>
void
exportDat(const TTree &tree, std::string outName)
{
    unsigned int i =0;
    double thickness = 1;
    std::ofstream fout;
    fout.open(outName.c_str());
    for (auto s : tree.myVectSegments) {
      // test if the segment is the root or its parent we do not display (already done).
      if (s.myIndex == 0 || s.myIndex == 1)
        continue;
      DGtal::Z3i::RealPoint distal = s.myCoordinate;
      DGtal::Z3i::RealPoint proxital = tree.myVectSegments[tree.myVectParent[s.myIndex]].myCoordinate;
      fout << distal[0] << " " << distal[1] << " " << distal[2] << " ";
      fout << proxital[0] << " " << proxital[1] << " " << proxital[2] << " ";
      fout << tree.myVectSegments[s.myIndex].myRadius*thickness << std::endl;
      i++;
    }
    fout.close();
}

#ifdef WITH_VISU3D_QGLVIEWER
template<typename TTree>
int
display3DTree(const TTree &tree )
{
    int arg = 0;
      QApplication application(arg,NULL);
      typedef DGtal::Viewer3D<> MyViewer;
      MyViewer viewer;
      viewer.show();
      unsigned int i = 0;
      double thickness = 1;
      viewer << DGtal::CustomColors3D(DGtal::Color(0,0,250),DGtal::Color(0,0,250));
      DGtal::Z3i::RealPoint p1 = tree.myVectSegments[1].myCoordinate;
      DGtal::Z3i::RealPoint p2 = tree.myVectSegments[tree.myVectParent[tree.myVectSegments[1].myIndex]].myCoordinate;
      viewer.addBall(p2,tree.myVectSegments[1].myRadius);
      viewer << DGtal::CustomColors3D(DGtal::Color(0,250,0),DGtal::Color(0,250,0));
      viewer.addCylinder(p1,p2,tree.myVectSegments[1].myRadius*thickness);
      
      for (auto s : tree.myVectSegments) {
        // test if the segment is the root or its parent we do not display (already done).
        if (s.myIndex == 0 || s.myIndex == 1)
          continue;
        DGtal::Z3i::RealPoint distal = s.myCoordinate;
        DGtal::Z3i::RealPoint proxital = tree.myVectSegments[tree.myVectParent[s.myIndex]].myCoordinate;
        viewer << DGtal::CustomColors3D(DGtal::Color(250,0,0),DGtal::Color(250,0,0));
        viewer.addBall(distal,tree.myVectSegments[s.myIndex].myRadius);
        viewer << DGtal::CustomColors3D(DGtal::Color(0,250,0),DGtal::Color(0,250,0));
        //viewer.addBall(distal,1);
        viewer.addCylinder(distal,proxital,tree.myVectSegments[s.myIndex].myRadius*thickness);
        i++;
      }
      /*
       //Display Sphere domaine
       viewer << DGtal::CustomColors3D(DGtal::Color(0,0,250,10),DGtal::Color(0,0,250,10));
       viewer.addBall(DGtal::Z3i::RealPoint(0,0,0),my_rPerf);
       */
      viewer<< MyViewer::updateDisplay;
      return application.exec();
    }
#endif


int main(int argc, char **argv)
{
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  app.description("Generated a 3D tree using the CCO algorithm. By default it generates a 3D mesh.");
  int nbTerm {500};
  double aPerf {20000};
  bool verbose {false};
  bool display3D {false};
  std::string nameImgDom {""};
  std::string outputMeshName {"result.off"};
  std::string exportDatName {""};
  std::string exportXMLName {""};
    
  std::vector<int> postInitV {-1,-1};
  app.add_option("-n,--nbTerm,1", nbTerm, "Set the number of terminal segments.", true);
  app.add_option("-a,--aPerf,2", aPerf, "The value of the input parameter A perfusion.", true);
  app.add_option("--organDomain,-d", nameImgDom, "Define the organ domain using a mask image (organ=255).");
  app.add_option("-o,--outputName", outputMeshName, "Output the 3D mesh", true);
  app.add_option("-e,--export", exportDatName, "Output the 3D mesh", true);
  app.add_option("-x,--exportXML", exportXMLName, "Output the resulting gaph as xml file", true);
  
#ifdef WITH_VISU3D_QGLVIEWER
  app.add_flag("--view", display3D, "display 3D view using QGLViewer");
#endif
  auto pInit = app.add_option("-p,--posInit", postInitV, "Initial position of root, if not given the position of point is determined from the image center")
  ->expected(3);
  app.add_flag("-v,--verbose", verbose);
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  DGtal::Z3i::Point ptRoot(postInitV[0], postInitV[1], postInitV[2]);

  if(nameImgDom != "" ){
    typedef ImageMaskDomainCtrl<3> TImgContrl;
    typedef  CoronaryArteryTree<TImgContrl, 3> TTree;
    TImgContrl aDomCtr (nameImgDom, 128, 100);
    TImgContrl::TPointI pM = aDomCtr.maxDistantPointFromBorder();
    unsigned int distSearchRoot =  static_cast< unsigned int >(aDomCtr.myDistanceImage(pM) /2.0);

    TTree tree  (aPerf, nbTerm, 1.0, pM);
    tree.myDomainController = aDomCtr;
    constructTreeMaskDomain(tree, distSearchRoot, verbose);
    
    XMLHelpers::writeTreeToXml(tree, "tree_3D.xml");
    exportResultingMesh(tree, outputMeshName);
    #ifdef WITH_VISU3D_QGLVIEWER
    if (display3D) display3DTree(tree);
    #endif
    if (exportDatName != "") exportDat(tree, exportDatName);
    if (exportXMLName != "") XMLHelpers::writeTreeToXml(tree, exportXMLName.c_str());
  }
  else
  {
    typedef CircularDomainCtrl<3> SphereDomCtrl;
    typedef  CoronaryArteryTree<SphereDomCtrl, 3> TTree;
    SphereDomCtrl::TPoint pCenter (0,0,0);
    TTree tree  (aPerf, nbTerm, 1.0, pCenter);
    constructTreeImplicitDomain(tree, verbose);
    XMLHelpers::writeTreeToXml(tree, "tree_3D.xml");
    exportResultingMesh(tree, outputMeshName);
    #ifdef WITH_VISU3D_QGLVIEWER
    if (display3D) display3DTree(tree);
    #endif
    if (exportDatName != "") exportDat(tree, exportDatName);
    if (exportXMLName != "") XMLHelpers::writeTreeToXml(tree,
                                                          exportXMLName.c_str());
  }
  return EXIT_SUCCESS;
}
