/**
 *  generateTree3D: main program to generate 2D tree from OpenCCO implementation
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
constructTreeMaskDomain(TTree &aTree,
                       bool verbose)
{
    clock_t start, end;
    start = clock();
    ExpandTreeHelpers::initFirstElemTree(aTree, verbose);
    ExpandTreeHelpers::expandTree(aTree, verbose);
    end = clock();
    printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
}


template<typename TTree>
void
constructTreeImplicitDomain(TTree &aTree, std::string outputMeshName,
                            std::string exportXMLName,
                            std::string exportDatName, bool verbose,
                            bool display3D)
{
    clock_t start, end;
    start = clock();
    ExpandTreeHelpers::initFirstElemTree(aTree, verbose);
    ExpandTreeHelpers::expandTree(aTree, verbose);
    end = clock();
    printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
    if (exportDatName != "") exportDat(aTree, exportDatName);
    if (exportXMLName != "") XMLHelpers::writeTreeToXml(aTree,
                                                        exportXMLName.c_str());

    XMLHelpers::writeTreeToXml(aTree, "tree_3D.xml");
    exportResultingMesh(aTree, outputMeshName);
    #ifdef WITH_VISU3D_QGLVIEWER
        if (display3D) display3DTree(aTree);
    #endif
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
      
      viewer<< MyViewer::updateDisplay;
      return application.exec();
    }
#endif


int main(int argc, char **argv)
{
  srand ((int) time(NULL));
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  app.description("Generated a 3D tree using the CCO algorithm. By default it generates a 3D mesh.");
  int nbTerm {500};
  double aPerf {20000};
  double gamma {3.0};

  double minDistanceToBorder {5.0};
  bool verbose {false};
  bool display3D {false};
    
  std::string nameImgDom {""};
  std::string outputMeshName {"result.off"};
  std::string exportDatName {""};
  std::string exportXMLName {""};
  std::vector<int> postInitV {-1,-1,-1};
  bool squaredImplDomain {false};

  app.add_option("-n,--nbTerm,1", nbTerm, "Set the number of terminal segments.", true);
  app.add_option("-a,--aPerf,2", aPerf, "The value of the input parameter A perfusion.", true);
  app.add_option("-g,--gamma", gamma, "The value of the gamma parameter.", true);
  app.add_option("--organDomain,-d", nameImgDom, "Define the organ domain using a mask image (organ=255).");
  app.add_option("-m,--minDistanceToBorder", minDistanceToBorder, "Set the minimal distance to border. Works only  with option organDomain else it has not effect", true);

  app.add_option("-o,--outputName", outputMeshName, "Output the 3D mesh", true);
  app.add_option("-e,--export", exportDatName, "Output the 3D mesh", true);
  app.add_option("-x,--exportXML", exportXMLName, "Output the resulting gaph as xml file", true);
  app.add_flag("-s,--squaredDom",squaredImplDomain , "Use a squared implicit domain instead a sphere (is used only without --organDomain)");

  auto pInit = app.add_option("-p,--posInit", postInitV, "Initial position of root, if not given the position of point is determined from the image center")
    ->expected(3);

#ifdef WITH_VISU3D_QGLVIEWER
  app.add_flag("--view", display3D, "display 3D view using QGLViewer");
#endif
  app.add_flag("-v,--verbose", verbose);
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  DGtal::Z3i::Point ptRoot(postInitV[0], postInitV[1], postInitV[2]);

  if(nameImgDom != "" ){
    typedef ImageMaskDomainCtrl<3> TImgContrl;
    typedef  CoronaryArteryTree<TImgContrl, 3> TTree;
    TImgContrl aDomCtr;
    PointI<3> pM;
      if (!pInit->empty())
      {
          pM[0] = postInitV[0];
          pM[1] = postInitV[1];
          pM[2] = postInitV[2];
          aDomCtr = TImgContrl(nameImgDom, 128, pM, 100);
      }
      else
      {
          aDomCtr = TImgContrl(nameImgDom, 128, 100);
      }
    
    aDomCtr.myMinDistanceToBorder = minDistanceToBorder;
    TTree tree  (aPerf, nbTerm, aDomCtr, 1.0);
    tree.my_gamma = gamma;

    constructTreeMaskDomain(tree, verbose);
    
    XMLHelpers::writeTreeToXml(tree, "tree_3D.xml");
    exportResultingMesh(tree, outputMeshName);
    #ifdef WITH_VISU3D_QGLVIEWER
    if (display3D) display3DTree(tree);
    #endif
    if (exportDatName != "") exportDat(tree, exportDatName);
    if (exportXMLName != "") XMLHelpers::writeTreeToXml(tree, exportXMLName.c_str());
  }
  else if (squaredImplDomain)
  {
      typedef SquareDomainCtrl<3> SqDomCtrl;
      typedef  CoronaryArteryTree<SqDomCtrl, 3> TTree;
      PointD<3> pCenter (0,0,0);
      SqDomCtrl aCtr(1.0 ,pCenter);
      TTree tree  (aPerf, nbTerm, aCtr,  1.0);
      tree.my_gamma = gamma;
      constructTreeImplicitDomain(tree, outputMeshName,
                                    exportXMLName,
                                    exportDatName, verbose, display3D);
  }
  else
  {
    typedef CircularDomainCtrl<3> SphereDomCtrl;
    typedef  CoronaryArteryTree<SphereDomCtrl, 3> TTree;
    PointD<3> pCenter (0,0,0);
    SphereDomCtrl aCtr(1.0 ,pCenter);
    TTree tree  (aPerf, nbTerm, aCtr,  1.0);
    tree.my_gamma = gamma;
    constructTreeImplicitDomain(tree, outputMeshName,
                                  exportXMLName,
                                  exportDatName, verbose, display3D);    
  }
  return EXIT_SUCCESS;
}
