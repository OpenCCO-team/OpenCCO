#include "CLI11.hpp"
#include "TreeImageRenderer.h"

#include "DGtal/io/colormaps/GradientColorMap.h"


void animation(const std::string & radii_filename, const std::string & vertices_filename, const std::string & edges_filename, int duration)
{
	// renderer initialized from files
	TreeImageRenderer<2> renderer(radii_filename, vertices_filename, edges_filename);

	renderer.treeConstructionAnimation("anim.svg", duration * 1000);
}

void test()
{
	TreeImageRenderer<2>::TPoint p1(0, 0);
	TreeImageRenderer<2>::TPoint p2(10, 10);
	TreeImageRenderer<2>::TDomain d(p1, p2);
	TreeImageRenderer<2>::TImage image(d);

	TreeImageRenderer<2>::TPointD line_start(7.8, 1.7);
	TreeImageRenderer<2>::TPointD line_end(2.99, 5.67);

	TreeImageRenderer<2>::TPoint line_start_r(7, 1);
	TreeImageRenderer<2>::TPoint line_end_r(7, 1);	

	drawBresenhamLine<2>(image, line_start_r, line_end_r);
}

int main(int argc, char *const *argv)
{
	// fixed dimension, for now we only animate in 2D
	const int dimension = 2;

	// parse command line using CLI ----------------------------------------------
	CLI::App app;
	std::string radii_filename = "radius.dat";
	std::string vertices_filename = "vertex.dat";
	std::string edges_filename = "edges.dat";
	double duration = 10.0;		// base duration of 10s
	std::string output_filename = "CCOanimation.svg";

	app.add_option("-r,--radii", radii_filename, "File containing the radii of the vertices.");
	app.add_option("-v,--vertices", vertices_filename, "File containing the coordinates of the vertices.");
	app.add_option("-e,--edges", edges_filename, "File containing the edges data.");
	app.add_option("-d,--duration", duration, "Duration of the animation, in milliseconds.")
		->check(CLI::PositiveNumber);	// duration needs to be positive
	app.add_option("-o,--output", output_filename, "File to write the animation to.");

	app.get_formatter()->column_width(40);
	CLI11_PARSE(app, argc, argv);
	// END parse command line using CLI ----------------------------------------------

	//animation(radii_filename, vertices_filename, edges_filename, duration);

	test();

	return 0;
}