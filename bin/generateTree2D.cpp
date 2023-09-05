/**
 *  generateTree2D: main program to generate 2D tree from OpenCCO implementation 
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

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "CLI11.hpp"
#include "XmlHelpers.h"

#include "CoronaryArteryTree.h"
#include "GeomHelpers.h"
#include "ExpandTreeHelpers.h"

#include "XmlHelpers.h"


/**
 * Function to construct the tree with the help ConstructionHelpers by using a domain of reconstruction defined fram an image (ImageMaskDomainCtrl)l
 */
template<typename TTree>
void
constructTreeMaskDomain(TTree &aTree, bool verbose)
{
	clock_t start, end;
	start = clock();
	ExpandTreeHelpers::initFirstElemTree(aTree);
	ExpandTreeHelpers::expandTree(aTree, verbose);
	end = clock();
	printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
}


template<typename TTree>
void
constructTreeImplicitDomain(TTree &aTree,std::string exportXMLName, bool verbose)
{
    clock_t start, end;
    start = clock();
    ExpandTreeHelpers::initFirtElemTree(aTree, verbose);
    ExpandTreeHelpers::expandTree(aTree);
    end = clock();
    printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
    XMLHelpers::writeTreeToXml(aTree, "tree_2D.xml");
    if (exportXMLName != "") XMLHelpers::writeTreeToXml(aTree,
                                                        exportXMLName.c_str());
    std::string filename = "testCCO_"+std::to_string(aTree.my_NTerm)+".eps";
    aTree.exportBoardDisplay(filename.c_str(), 1.0);
    aTree.exportBoardDisplay("result.eps", 1.0);
    aTree.exportBoardDisplay("result.svg", 1.0);
    aTree.myBoard.clear();
}



/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{
	clock_t start, end;
	srand ((int) time(NULL));

	// parse command line using CLI ----------------------------------------------
	CLI::App app;
	int nbTerm {1000};
	double aPerf {20000};
	bool verbose {false};
	bool squaredImplDomain {false};
	double gamma {3.0};
	double minDistanceToBorder {5.0};
	std::string nameImgDom {""};
	std::vector<int> postInitV {-1,-1};
	std::string exportXMLName {""};
	
	app.add_option("-n,--nbTerm,1", nbTerm, "Set the number of terminal segments.", true);
	app.add_option("-a,--aPerf,2", aPerf, "The value of the input parameter A perfusion.", true);
	app.add_option("-g,--gamma", gamma, "The value of the gamma parameter.", true);
	app.add_option("-m,--minDistanceToBorder", minDistanceToBorder, "Set the minimal distance to border. Works only  with option organDomain else it has not effect", true);
	app.add_option("--organDomain,-d", nameImgDom, "Define the organ domain using a mask image (organ=255).");
	app.add_option("-x,--exportXML", exportXMLName, "Output the resulting gaph as xml file", true);
	app.add_flag("-s,--squaredDom",squaredImplDomain , "Use a squared implicit domain instead a sphere (is used only without --organDomain)");
	auto pInit = app.add_option("-p,--posInit", postInitV, "Initial position of root, if not given the position of point is determined from the image center")
		->expected(2);
	app.add_flag("-v,--verbose", verbose);
	app.get_formatter()->column_width(40);
	CLI11_PARSE(app, argc, argv);
	// END parse command line using CLI ----------------------------------------------
	
	DGtal::Z2i::Point ptRoot(postInitV[0], postInitV[1]);
  
	if(nameImgDom != ""){
		start = clock();
		typedef ImageMaskDomainCtrl<2> TImgContrl;
		typedef  CoronaryArteryTree<TImgContrl, 2> TTree;
		TImgContrl aDomCtr;
		PointI<2> pM;
		if (!pInit->empty())
		{
			pM[0] = postInitV[0];
			pM[1] = postInitV[1];
			aDomCtr = TImgContrl(nameImgDom, 128, pM, 100);
		}
		else
		{
			aDomCtr = TImgContrl(nameImgDom, 128, 100);
		}
    // parse command line using CLI ----------------------------------------------
    CLI::App app;
    int nbTerm {1000};
    double aPerf {20000};
    bool verbose {false};
    bool squaredImplDomain {false};
    double gamma {3.0};
    double minDistanceToBorder {5.0};
    std::string nameImgDom {""};
    std::vector<int> postInitV {-1,-1};
    std::string exportXMLName {""};
    std::string outputNameEPS {"result.eps"};
    std::string outputNameSVG {"result.svg"};
    
    app.add_option("-n,--nbTerm,1", nbTerm, "Set the number of terminal segments.", true);
    app.add_option("-a,--aPerf,2", aPerf, "The value of perfusion area.", true);
    app.add_option("-g,--gamma", gamma, "The value of the gamma parameter.", true);
    app.add_option("-m,--minDistanceToBorder", minDistanceToBorder, "Set the minimal distance to border. Works only  with option organDomain else it has not effect", true);
    app.add_option("-o,--outputEPS", outputNameEPS, "Output the result into EPS format", true);
    app.add_option("-e,--exportSVG", outputNameSVG, "Export the result into SVG format", true);
    app.add_option("--organDomain,-d", nameImgDom, "Define the organ domain using a mask image (organ=255).");
    app.add_option("-x,--exportXML", exportXMLName, "Output the resulting graph as xml file", true);
    app.add_flag("-s,--squaredDom",squaredImplDomain , "Use a squared implicit domain instead a sphere (is used only without --organDomain)");
    auto pInit = app.add_option("-p,--posInit", postInitV, "Initial position of root, if not given the position of point is determined from the image center")
    ->expected(2);
    app.add_flag("-v,--verbose", verbose);
    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------------
    
    DGtal::Z2i::Point ptRoot(postInitV[0], postInitV[1]);
  
    if(nameImgDom != ""){
        start = clock();
        typedef ImageMaskDomainCtrl<2> TImgContrl;
        typedef  CoronaryArteryTree<TImgContrl, 2> TTree;
        TImgContrl aDomCtr;
        TImgContrl::TPointI pM;
        if (!pInit->empty())
        {
            pM[0] = postInitV[0];
            pM[1] = postInitV[1];
            aDomCtr = TImgContrl(nameImgDom, 128, pM, 100);
        }
        else
        {
            aDomCtr = TImgContrl(nameImgDom, 128, 100);
        }
        aDomCtr.myMinDistanceToBorder = minDistanceToBorder;
        TTree tree  (aPerf, nbTerm, aDomCtr,  1.0);
        tree.my_gamma = gamma;
        constructTreeMaskDomain(tree, verbose);
        XMLHelpers::writeTreeToXml(tree, "tree_2D.xml");
        if (exportXMLName != "") XMLHelpers::writeTreeToXml(tree, exportXMLName.c_str());
        
        std::string filename = "testCCO_"+std::to_string(nbTerm)+".eps";
        tree.exportBoardDisplay(filename.c_str(), 1.0);
        tree.exportBoardDisplay(outputNameEPS.c_str(), 1.0);
        tree.exportBoardDisplay(outputNameSVG.c_str(), 1.0);
        tree.myBoard.clear();
        end = clock();
        printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
        
    }
    else if (squaredImplDomain)
    {
        start = clock();
        typedef SquareDomainCtrl<2> SqDomCtrl;
        typedef  CoronaryArteryTree<SqDomCtrl, 2> TTree;
        SqDomCtrl::TPoint pCenter (0,0);
        SqDomCtrl aCtr(1.0,pCenter);
        TTree tree  (aPerf, nbTerm, aCtr, 1.0);
        tree.my_gamma = gamma;
        constructTreeImplicitDomain(tree, exportXMLName, verbose);
        end = clock();
        printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
    }else {
        start = clock();
        typedef CircularDomainCtrl<2> DiskDomCtrl;
        typedef  CoronaryArteryTree<DiskDomCtrl, 2> TTree;
        DiskDomCtrl::TPoint pCenter (0,0);
        DiskDomCtrl aCtr(1.0,pCenter);
        TTree tree  (aPerf, nbTerm, aCtr, 1.0);
        tree.my_gamma = gamma;
        constructTreeImplicitDomain(tree, exportXMLName, verbose);
        end = clock();
        printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);

		aDomCtr.myMinDistanceToBorder = minDistanceToBorder;
		TTree tree  (aPerf, nbTerm, aDomCtr,  1.0);
		tree.my_gamma = gamma;
		constructTreeMaskDomain(tree, verbose);
		XMLHelpers::writeTreeToXml(tree, "tree_2D.xml");
		if (exportXMLName != "") XMLHelpers::writeTreeToXml(tree, exportXMLName.c_str());
		
		std::string filename = "testCCO_"+std::to_string(nbTerm)+".eps";
		tree.exportBoardDisplay(filename.c_str(), 1.0);
		tree.exportBoardDisplay("result.eps", 1.0);
		tree.exportBoardDisplay("result.svg", 1.0);
		tree.myBoard.clear();
		end = clock();
		printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
		
	}
	else if (squaredImplDomain)
	{
		start = clock();
		typedef SquareDomainCtrl<2> SqDomCtrl;
		typedef  CoronaryArteryTree<SqDomCtrl, 2> TTree;
		PointD<2> pCenter (0,0);
		SqDomCtrl aCtr(1.0,pCenter);
		TTree tree  (aPerf, nbTerm, aCtr, 1.0);
		tree.my_gamma = gamma;
		constructTreeImplicitDomain(tree, exportXMLName, verbose);
		end = clock();
		printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
	}else {
		start = clock();
		typedef CircularDomainCtrl<2> DiskDomCtrl;
		typedef  CoronaryArteryTree<DiskDomCtrl, 2> TTree;
		PointD<2> pCenter (0,0);
		DiskDomCtrl aCtr(1.0,pCenter);
		TTree tree  (aPerf, nbTerm, aCtr, 1.0);
		tree.my_gamma = gamma;
		constructTreeImplicitDomain(tree, exportXMLName, verbose);
		end = clock();
		printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);

	}
	
	return EXIT_SUCCESS;
}
