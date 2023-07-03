#include "CLI11.hpp"
#include "TreeImageRenderer.h"


#include <iomanip>
#include <iostream>
#include <random>
#include <chrono>
#include <map>

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

void test()
{
	// donut renderer :)
	
	int scale_factor = 60;

	// tore radii
	double r = 1.0 * scale_factor;
	double R = 2 * r;

	TPoint<3> p1(-4 * r, -4 * r, -2 * r);
	TPoint<3> p2(-1 * p1);

	TImage<3> torvol(TDomain<3>(p1, p2));

	for(const TPoint<3> & p : torvol.domain())
	{
		double var = sqrt(p[0] * p[0] + p[1] * p[1]) - R;
		double eq_tore_p = var * var + p[2] * p[2] - r * r;

		if (eq_tore_p <= 0.0)
		{
			// p belongs to the tore volume
			torvol.setValue(p, 255);
		}
		else
		{
			torvol.setValue(p, 0.0);
		}
	}

	saveRender<3>(torvol, "tore");
}

int main(int argc, char *const *argv)
{
	// parse command line using CLI ----------------------------------------------
	CLI::App app;
	unsigned int output_width = 0;
	std::string domain_filename = "";
	std::string radii_filename = "radius.dat";
	std::string vertices_filename = "vertex.dat";
	std::string edges_filename = "edges.dat";
	std::string output_filename = "realistic_render";
	double sigma = 3.0;

	auto dom_group = app.add_option_group("Render domain");
	dom_group->add_option("-w,--width", output_width, "Width of the output image, in pixels. Aspect ratio is constrained by the position of the points.", true)
		->check(CLI::Range(50, 10000));			// nobody would create an image with a width outside this range, right ?
	dom_group->add_option("-d,--domain", domain_filename, "Image file defining the organ domain (organ >= 128).");
	dom_group->require_option(1); 				// mandatory to use one (not less, not more) of these two options

	app.add_option("-r,--radii", radii_filename, "File containing the radii of the vertices.");
	app.add_option("-v,--vertices", vertices_filename, "File containing the coordinates of the vertices.");
	app.add_option("-e,--edges", edges_filename, "File containing the edges data.");
	app.add_option("-o,--output", output_filename, "File to write the animation to (without file extension).");
	app.add_option("-s,--sigma", sigma, "The standard deviation of the noise for the realistic render.")
		->check(CLI::PositiveNumber);			// SNR is strictly positive

	app.get_formatter()->column_width(40);
	CLI11_PARSE(app, argc, argv);
	// END parse command line using CLI ----------------------------------------------

	
	int dimension = checkDimension(vertices_filename);


	if (dimension == 2)
	{
		TImage<2> img{ TDomain<2>() };

		// renderer initialized from files
		if(domain_filename == "")
		{
			TreeImageRenderer<2> renderer(radii_filename, vertices_filename, edges_filename);

			// output_width is ignored if domain_filename was defined
			img = renderer.realisticRender(sigma, output_width);
		}
		else
		{
			TreeImageRenderer<2> renderer(radii_filename, vertices_filename, edges_filename, domain_filename);

			// output_width is ignored if domain_filename was defined
			img = renderer.realisticRender(sigma, output_width);
		}

		saveRender<2>(img, output_filename);
	}
	else if (dimension == 3)
	{
		TImage<3> img{ TDomain<3>() };

		// renderer initialized from files
		if(domain_filename == "")
		{
			TreeImageRenderer<3> renderer(radii_filename, vertices_filename, edges_filename);

			// output_width is ignored if domain_filename was defined
			img = renderer.realisticRender(sigma, output_width);
		}
		else
		{
			std::cout << "here" << std::endl;
			TreeImageRenderer<3> renderer(radii_filename, vertices_filename, edges_filename, domain_filename);

			// output_width is ignored if domain_filename was defined
			img = renderer.realisticRender(sigma, output_width);
		}

		saveRender<3>(img, output_filename);
	}

	//test();

	return 0;
}