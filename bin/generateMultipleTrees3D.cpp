#include <iostream>

#include "CLI11.hpp"

#include "CohabitingTrees.h"
#include "ExpandTreeHelpers.h"


template<int TDim>
void templa()
{
	CircularDomainCtrl<TDim> circ_domctr(1.0, PointD<TDim>::diagonal(0.0));

	// data vectors
	double perfs = 10000;
	std::vector<unsigned int> terms {50, 50};

	CohabitingTrees<CircularDomainCtrl<TDim>, TDim> ctree(perfs, terms, circ_domctr);

	ExpandTreeHelpers::expandCohabitingTrees(ctree, 200);

	ctree.writeTreesToXML();
}

void templamaks()
{
	std::string domain_file = "Samples/shape.png";
	CircularDomainCtrl<2> circ_domctr(1.0, PointD<2>::diagonal(0.0));

	// data vectors
	double perfs = 100;
	std::vector<unsigned int> terms {4, 4};

	CohabitingTrees<CircularDomainCtrl<2>, 2> ctree(perfs, terms, circ_domctr);

	ExpandTreeHelpers::expandCohabitingTrees(ctree, 200);

	ctree.exportTreesDisplays();
}

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
	app.add_option("-a,--aPerf,2", aPerf, "The value of perfusion volume.", true);
	app.add_option("-g,--gamma", gamma, "The value of the gamma parameter.", true);
	app.add_option("--organDomain,-d", nameImgDom, "Define the organ domain using a mask image (organ=255).");
	app.add_option("-m,--minDistanceToBorder", minDistanceToBorder, "Set the minimal distance to border. Works only with option --organDomain else it has no effect", true);
	app.add_option("-o,--outputName", outputMeshName, "Output the 3D mesh into OFF format", true);
	app.add_option("-e,--export", exportDatName, "Output the 3D mesh into text file", true);
	app.add_option("-x,--exportXML", exportXMLName, "Output the resulting graph as xml file", true);
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

	return 0;
}