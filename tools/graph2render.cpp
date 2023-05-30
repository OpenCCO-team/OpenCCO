#include "TreeImageRenderer.h"

#include "DGtal/io/writers/STBWriter.h"
#include <DGtal/io/colormaps/GradientColorMap.h>

#include "CLI11.hpp"

/*
void createPerlinNoiseBackground(unsigned int width, unsigned int height)
{
	TPoint p1(0, 0);
	TPoint p2(width, height);
	TDomain domain(p1, p2);
	TImage2D background(domain);

	DGtal::HueShadeColorMap<double>  hueColorMap (-50.0, 50.0);

	DGtal::STBWriter<TImage2D, DGtal::HueShadeColorMap<double>>::exportPNG("aFilename.png", background, hueColorMap);

	return;
}
*/

void saveImage(const TImage2D & image, const std::string & filename)
{
	auto min_val = std::min_element(image.constRange().begin(), image.constRange().end());
	auto max_val = std::max_element(image.constRange().begin(), image.constRange().end());

	DGtal::GradientColorMap<double> gradient_cmap(*min_val, *max_val);
	gradient_cmap.addColor(DGtal::Color::Black);
    gradient_cmap.addColor(DGtal::Color::White);


    DGtal::STBWriter<TImage2D, DGtal::GradientColorMap<double>>::exportPNG(filename, image, gradient_cmap);
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

    TreeImageRenderer renderer(output_width, radii_filename, vertices_filename, edges_filename);

    renderer.createTreeImage();

    renderer.createDistanceMap();

    saveImage(renderer.treeImage(), "treeimage.png");
    saveImage(renderer.distanceMap(), "distancemap.png");
}