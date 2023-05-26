#include "TreeImageRenderer.h"
#include "DGtal/io/writers/GenericWriter.h"

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

int main(int argc, char *const *argv)
{
	// parse command line using CLI ----------------------------------------------
    CLI::App app;
    unsigned int output_width{1000};

    app.add_option("-w,--width", output_width, "Width of the output image, in pixels. Aspect ratio is constrained by the position of the points.", true);

    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------------

    TreeImageRenderer renderer(output_width);

    //renderer.createDistanceAlphaChannel();
}