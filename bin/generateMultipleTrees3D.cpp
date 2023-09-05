#include <iostream>
#include <tuple>
#include <chrono>

#include "CLI11.hpp"

#include "CohabitingTrees.h"
#include "ExpandTreeHelpers.h"


template<typename CohabTrees>
void constructTrees(CohabTrees & trees, std::string outname)
{
	const auto start = std::chrono::steady_clock::now();
	ExpandTreeHelpers::expandCohabitingTrees(trees);
	const auto end = std::chrono::steady_clock::now();

	const std::chrono::duration<double> expansion_time = end - start;
	std::cout << "Execution time: " << expansion_time.count() << " sec" << std::endl;

	trees.writeTreesToXML(outname);
}


int main(int argc, char **argv)
{
	srand ((int) time(NULL));

	// parse command line using CLI ----------------------------------------------
	CLI::App app;
	app.description("Generates multiple 3D tree using the CCO algorithm.");
	int nbTree {2};
	int nbTerm {500};
	double aPerf {20000};
	double gamma {3.0};

	double minDistanceToBorder {5.0};
	bool verbose {false};
	bool display3D {false};
	std::string nameImgDom {""};
	std::string outputMeshName {"result.off"};
	std::string exportDatName {""};
	std::string exportXMLName {"tree"};
	std::vector< std::tuple<int, int , int> > posInitV;
	bool squaredImplDomain {false};

	app.add_option("-t,--nbTree", nbTree, "Set the number of vascular trees.", true);
	app.add_option("-n,--nbTerm,1", nbTerm, "Set the number of terminal segments.", true);
	app.add_option("-a,--aPerf,2", aPerf, "The value of perfusion volume.", true);
	app.add_option("-g,--gamma", gamma, "The value of the gamma parameter.", true);
	app.add_option("-d,--organDomain", nameImgDom, "Define the organ domain using a mask image (organ=255)");
	app.add_option("-m,--minDistanceToBorder", minDistanceToBorder, "Set the minimal distance to border. Works only with option --organDomain else it has no effect", true);
	app.add_option("-o,--outputName", outputMeshName, "Output the 3D mesh into OFF format", true);
	app.add_option("-e,--export", exportDatName, "Output the 3D mesh into text file", true);
	app.add_option("-x,--exportXML", exportXMLName, "Output the resulting graph as xml file", true);
	app.add_flag("-s,--squaredDom",squaredImplDomain , "Use a squared implicit domain instead a sphere (is used only without --organDomain)");
	auto pInit = app.add_option("-p,--posInit", posInitV, "Initial position of roots (x,y,z ... x,y,z etc), if not given the positions are evenly spread in the organ domain");

	#ifdef WITH_VISU3D_QGLVIEWER
	app.add_flag("--view", display3D, "display 3D view using QGLViewer");
	#endif
	app.add_flag("-v,--verbose", verbose);
	app.get_formatter()->column_width(40);
	CLI11_PARSE(app, argc, argv);
	// END parse command line using CLI ----------------------------------------------

	if(!posInitV.empty() && (posInitV.size() != nbTree))
	{
		std::cout << "Input error : The number of initial positions of roots must be equal to the number of vascular trees or be unspecifed." << std::endl;
		return -1;
	}

	std::vector< PointI<3> > starting_points;
	for(int i = 0; i < posInitV.size(); i++)
	{
		starting_points.emplace_back(std::get<0>(posInitV[i]),
									 std::get<1>(posInitV[i]),
									 std::get<2>(posInitV[i]));
	}

	if(nameImgDom.empty())
	{
		if(squaredImplDomain)	// square implicit domain
		{
			SquareDomainCtrl<3> square_domain(1.0, PointD<3>(0,0,0));
			std::vector<unsigned int> terms(nbTree, nbTerm);

			CohabitingTrees<SquareDomainCtrl<3>, 3> ctree(aPerf, terms, gamma, starting_points, square_domain);

			constructTrees(ctree, exportXMLName);
		}
		else					// circular implicit domain
		{
			CircularDomainCtrl<3> circle_domain(1.0, PointD<3>(0,0,0));
			std::vector<unsigned int> terms(nbTree, nbTerm);

			CohabitingTrees<CircularDomainCtrl<3>, 3> ctree(aPerf, terms, gamma, starting_points, circle_domain);

			constructTrees(ctree, exportXMLName);
		}
	}
	else						// image mask domain
	{
		ImageMaskDomainCtrl<3> mask_domain(nameImgDom, 128);
		mask_domain.myMinDistanceToBorder = minDistanceToBorder;
		std::vector<unsigned int> terms(nbTree, nbTerm);

		CohabitingTrees<ImageMaskDomainCtrl<3>, 3> ctree(aPerf, terms, gamma, starting_points, mask_domain);
	
		constructTrees(ctree, exportXMLName);
	}

	return EXIT_SUCCESS;
}