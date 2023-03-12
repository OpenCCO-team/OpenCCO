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
    ExpandTreeHelpers::initFirtElemTree(aTree);
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
    ExpandTreeHelpers::initFirtElemTree(aTree, verbose);
    ExpandTreeHelpers::expandTree(aTree);
    end = clock();
    printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
}



/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{
    clock_t start, end;
    
    // parse command line using CLI ----------------------------------------------
    CLI::App app;
    int nbTerm {1000};
    double aPerf {20000};
    bool verbose {false};
    std::string nameImgDom {""};
    std::vector<int> postInitV {-1,-1};
    std::string exportXMLName {""};
    
    app.add_option("-n,--nbTerm,1", nbTerm, "Set the number of terminal segments.", true);
    app.add_option("-a,--aPerf,2", aPerf, "The value of the input parameter A perfusion.", true);
    app.add_option("--organDomain,-d", nameImgDom, "Define the organ domain using a mask image (organ=255).");
    app.add_option("-x,--exportXML", exportXMLName, "Output the resulting gaph as xml file", true);
    
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
        
        TTree tree  (aPerf, nbTerm, 1.0,aDomCtr.myCenter);
        tree.myDomainController = aDomCtr;
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
    else {
        start = clock();
        typedef CircularDomainCtrl<2> DiskDomCtrl;
        typedef  CoronaryArteryTree<DiskDomCtrl, 2> TTree;
        DiskDomCtrl::TPoint pCenter (0,0);
        TTree tree  (aPerf, nbTerm, 1.0, pCenter);
        DiskDomCtrl aCtr(tree.bParam.my_rPerf,pCenter);
        tree.myDomainController = aCtr;
        constructTreeImplicitDomain(tree, verbose);
        XMLHelpers::writeTreeToXml(tree, "tree_2D.xml");
        if (exportXMLName != "") XMLHelpers::writeTreeToXml(tree,
                                                            exportXMLName.c_str());
        std::string filename = "testCCO_"+std::to_string(nbTerm)+".eps";
        tree.exportBoardDisplay(filename.c_str(), 1.0);
        tree.exportBoardDisplay("result.eps", 1.0);
        tree.exportBoardDisplay("result.svg", 1.0);
        tree.myBoard.clear();
        end = clock();
        printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
        
    }
    
    return EXIT_SUCCESS;
}
