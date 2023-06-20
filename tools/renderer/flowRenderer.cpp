#include "TreeImageRenderer.h"
#include "CLI11.hpp"



/**
 * Returns the dimension of the space defined by the vertices.
 * @param vertices_filename the file containing the vertices coordinates.
 * @return The dimension.
 **/
int checkDimension(const std::string & vertices_filename)
{
	int dim = 0;

	std::ifstream v_file;

	v_file.open(vertices_filename);

	if(v_file.is_open())
	{
		// the first line contains the coordinates of the first point
		std::string line;
		getline(v_file, line);

		// separate the coordinates
		std::stringstream sline(line);

		std::istream_iterator<std::string> begin(sline);
		std::istream_iterator<std::string> end;

		std::vector<std::string> v_strings(begin, end);

		// amount of coordinates == dimension
		dim = v_strings.size();

		v_file.close();
	}

	return dim;
}




int main(int argc, char *const *argv)
{
	// parse command line using CLI ----------------------------------------------
	CLI::App app;
	unsigned int output_width = 1000;
	std::string radii_filename = "radius.dat";
	std::string vertices_filename = "vertex.dat";
	std::string edges_filename = "edges.dat";
	std::string output_filename = "flow_render";
	bool log_scale = false;

	app.add_option("-w,--width", output_width, "Width of the output image, in pixels. Aspect ratio is constrained by the position of the points.", true)
		->check(CLI::Range(100,10000));			// nobody would create images outside this range, right ?
	app.add_option("-r,--radii", radii_filename, "File containing the radii of the vertices.");
	app.add_option("-v,--vertices", vertices_filename, "File containing the coordinates of the vertices.");
	app.add_option("-e,--edges", edges_filename, "File containing the edges data.");
	app.add_option("-o,--output", output_filename, "File to write the animation to (format is relevant).");
	app.add_flag("-l,--log-scale", log_scale, "Uses a logarithmic scale of the flow instead of linear");

	app.get_formatter()->column_width(40);
	CLI11_PARSE(app, argc, argv);
	// END parse command line using CLI ----------------------------------------------

	int dimension = checkDimension(vertices_filename);

	if (dimension == 2)
	{
		// renderer initialized from files
		TreeImageRenderer<2> renderer(radii_filename, vertices_filename, edges_filename);

		TImage<2> img = renderer.flowRender(output_width);

		if(log_scale)
		{
			for(auto it = img.begin(); it != img.end(); it++)
			{
				*it = std::log1p(*it);
			}
		}

		saveRender<2>(img, output_filename);
	}
	else if (dimension == 3)
	{
		// renderer initialized from files
		TreeImageRenderer<3> renderer(radii_filename, vertices_filename, edges_filename);

		TImage<3> img = renderer.flowRender(output_width);

		if(log_scale)
		{
			for(auto it = img.begin(); it != img.end(); it++)
			{
				*it = std::log1p(*it);
			}
		}

		saveRender<3>(img, output_filename);
	}


	return 0;
}