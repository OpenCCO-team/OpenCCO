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

    app.add_option("-w,--width", output_width, "Width of the output image, in pixels. Aspect ratio is constrained by the position of the points.", true);
    app.add_option("-r,--radii", radii_filename, "File containing the radii of the vertices.");
    app.add_option("-v,--vertices", vertices_filename, "File containing the coordinates of the vertices.");
    app.add_option("-e,--edges", edges_filename, "File containing the edges data.");

    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------------

    int dimension = checkDimension(vertices_filename);

    if (dimension == 2)
    {
    	TreeImageRenderer<2> renderer(output_width, radii_filename, vertices_filename, edges_filename);

	    renderer.createTreeImage();

	    //renderer.createDistanceMap();

	    renderer.saveRender("render");
    }
    else if (dimension == 3)
    {
    	TreeImageRenderer<3> renderer(output_width, radii_filename, vertices_filename, edges_filename);

	    renderer.createTreeImage();

	    //renderer.createDistanceMap();

	    renderer.saveRender("render");
    }
    

    return 0;
}