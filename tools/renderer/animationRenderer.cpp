#include "CLI11.hpp"
#include "TreeImageRenderer.h"

#include "DGtal/io/colormaps/GradientColorMap.h"

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

	// renderer initialized from files
	TreeImageRenderer<2> renderer(radii_filename, vertices_filename, edges_filename);

	renderer.animationRender(output_filename, duration * 1000);

	return 0;
}