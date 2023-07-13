#include "TreeImageRenderer.h"

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/PGMWriter.h"
#include "DGtal/io/writers/STBWriter.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/io/writers/ITKWriter.h"

#include "GeomHelpers.h"

#include <random>
#include <limits>

#include "svg_elements.h"



///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////               OrganDomain<TDim>               ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


template<>
TreeImageRenderer<2>::OrganDomain::OrganDomain(const std::string & domain_filename)
 : myDomainMask(Domain<TSpace>()), isDefined(true)
{
	// if 2D image file when TDim = 3 : the image is in the first 2 dims, anything along z > 0 is zero
	// if 3D image file when TDim = 2 : the image is a slice of the 3D vol
	myDomainMask = DGtal::ITKReader< Image<2> >::importITK(domain_filename);
}



template<>
TreeImageRenderer<3>::OrganDomain::OrganDomain(const std::string & domain_filename)
 : myDomainMask(Domain<TSpace>()), isDefined(true)
{
	// if 2D image file when TDim = 3 : the image is in the first 2 dims, anything along z > 0 is zero
	// if 3D image file when TDim = 2 : the image is a slice of the 3D vol
	std::size_t pos = domain_filename.find_last_of('.');

	if(pos != std::string::npos)
	{
		const std::string ext = domain_filename.substr(pos+1);

		std::cout << "EXTENSION " << ext << " " << (ext == "vol") << std::endl;

		// vol is not supported by ITK, we handle it separately
		if(ext == "vol")
		{
			myDomainMask = DGtal::VolReader< Image<3> >::importVol(domain_filename);

			return;
		}
	}

	myDomainMask = DGtal::ITKReader< Image<3> >::importITK(domain_filename);
}


///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////          TreeImageRenderer<TDim>              ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// Constructor

template<int TDim>
TreeImageRenderer<TDim>::TreeImageRenderer(const std::string & radii_filename,
									 const std::string & vertices_filename,
									 const std::string & edges_filename,
									 const std::string & organ_dom_filename)
	: myOrganDomain(organ_dom_filename)
{
	importTreeData(radii_filename, vertices_filename, edges_filename);
}



template<int TDim>
TreeImageRenderer<TDim>::TreeImageRenderer(const std::string & radii_filename,
									 const std::string & vertices_filename,
									 const std::string & edges_filename)
	: myOrganDomain()
{
	importTreeData(radii_filename, vertices_filename, edges_filename);
}



// Methods

template<int TDim>
void TreeImageRenderer<TDim>::importTreeData(const std::string & radii_filename,
									   const std::string & vertices_filename,
									   const std::string & edges_filename)
{
	std::ifstream file;

	// get the vertex from vertex.dat
	file.open(vertices_filename);

	if(file.is_open())
	{
		std::string line;
		while(getline(file,line))
		{
			PointD<TDim> p;

			std::istringstream iss(line);

			for(unsigned int i = 0; i < TDim; i++)
			{
				iss >> p[i];
			}
			
			myTree.myPoints.push_back(p);
		}

		file.close();
	}
	else
	{
		std::cout << "Couldn't open " << vertices_filename << "." << std::endl;
	}

	// get the edge data from edges.dat
	file.open(edges_filename);

	if(file.is_open())
	{
		std::string line;
		
		// dump the first 3 lines
		getline(file,line);
		getline(file,line);
		getline(file,line);

		while(getline(file,line))
		{
			Segment s;

			std::istringstream iss(line);
			iss >> s.myProxitalIndex >> s.myDistalIndex >> s.myFlow >> s.myResistance;

			myTree.mySegments.push_back(s);
		}

		file.close();
	}
	else
	{
		std::cout << "Couldn't open " << edges_filename << "." << std::endl;
	}

	// get the radii from from radius.dat
	file.open(radii_filename);

	if(file.is_open())
	{
		std::string line;

		while(getline(file,line))
		{
			double r;
			std::istringstream iss(line);
			iss >> r;
			myTree.myRadii.push_back(r);
		}

		file.close();
	}
	else
	{
		std::cout << "Couldn't open " << radii_filename << "." << std::endl;;
	}
}



template<int TDim>
Image<TDim> TreeImageRenderer<TDim>::flowRender(unsigned int width)
{
	// domain
	Image<TDim> organ_img(createDomainImage(width, width/20));
	Image<TDim> flow_render(organ_img.domain());		// image of the same size

	// compute flow for each pixel
	for(const Segment & s : myTree.mySegments)
	{
		// variables useful for computation in subloop
		double dist_proxitaldistal = (myTree.myPoints[s.myDistalIndex] - myTree.myPoints[s.myProxitalIndex]).norm();
		double value = 7 + std::log1p(s.myFlow);

		// points bounding the segment
		std::vector< PointD<TDim> > v;
		for(int i = -1; i <= 1; i+=2)
		{
			v.emplace_back(myTree.myPoints[s.myProxitalIndex] 
							+ i * PointI<TDim>::diagonal(myTree.myRadii[s.myProxitalIndex]));
			v.emplace_back(myTree.myPoints[s.myDistalIndex]
							+ i * PointI<TDim>::diagonal(myTree.myRadii[s.myDistalIndex]));
		}

		// segment bounding box
		PointD<TDim> ub, lb;
		compBB<TDim>(v, ub, lb);

		for(const PointI<TDim> & p : Domain<TSpace>(lb, ub))
		{
			PointD<TDim> proj;
			bool isproj = projectOnStraightLine<TDim>(myTree.myPoints[s.myDistalIndex],
												myTree.myPoints[s.myProxitalIndex],
												p,
												proj);

			if(isproj)		// the projection belongs to the segment
			{
				double dist_pproj = (proj - p).norm();

				double dist_proxitalproj = (proj - myTree.myPoints[s.myProxitalIndex]).norm();
				double interpolated_radius = (myTree.myRadii[s.myDistalIndex] - myTree.myRadii[s.myProxitalIndex]) * dist_proxitalproj/dist_proxitaldistal
											+ myTree.myRadii[s.myProxitalIndex];

				if(dist_pproj < interpolated_radius && s.myFlow > flow_render(p))	// inside the segment and bigger flow	
				{
					flow_render.setValue(p, value);
				}
			}
			else		// the projection doesn't belong to the segment
			{			// the distance from p to the segment is either the distance to the distal or the distance to the proxital
						// offset by their respective radius
				double dist_pproxital = (myTree.myPoints[s.myProxitalIndex] - p).norm() - myTree.myRadii[s.myProxitalIndex];
				double dist_pdistal = (myTree.myPoints[s.myDistalIndex] - p).norm() - myTree.myRadii[s.myDistalIndex];

				double min_dist = std::min(dist_pproxital, dist_pdistal);

				if(min_dist < 0 && s.myFlow > flow_render(p))		// inside the segment and bigger flow
				{
					flow_render.setValue(p, value);
				}
			}
		}
	}

	normalizeImageValues<TDim>(flow_render, 0.0, 255.0);

	if(myOrganDomain.isDefined)
	{
		// blend the organ and the flow
		auto customblend = [](double a, double b) { return (b <= 0.001 ? a : b); };
		imageBlend<TDim>(organ_img, flow_render, customblend, flow_render);
	}

	return flow_render;
}



template<>
void TreeImageRenderer<2>::animationRender(const std::string & filename, int duration)
{
	// Right now the following code heavily utilizes the fact that the segment vector is sorted in a specific way.
	// If the logic of the CCO algorithm that produces the data used change the order of the segment for a reason or another
	// the following code might not output an animatoin that reproduces the logic of the CCO algorithm.
	//
	// Working on a way to preprocess the segment vector so that a good animation is guaranteed.	

	// copy segment data, they are modified by the animation process
	std::vector<Segment> segments = myTree.mySegments;

	// animation global variable
	int n = segments.size() / 2;						// number of terminal segments
	int segment_growth_dur = duration / n;				// duration, in milliseconds
	int delay = segment_growth_dur/4;					// intermediate duration for 2-part animations

	// find the min and the max of the flows
	double q_min = segments[0].myFlow;
	double q_max = q_min;
	for(const Segment & s : segments)
	{
		if(s.myFlow < q_min)
		{
			q_min = s.myFlow;
		}
		else if (s.myFlow > q_max)
		{
			q_max = s.myFlow;
		}
	}

	// color gradient for the vessels, depending on the value of the log of the flow
	DGtal::GradientColorMap<float> cmap_grad(std::log(q_min), std::log(q_max));
	cmap_grad.addColor( DGtal::Color( 79, 162, 198 ) );   // blue
	cmap_grad.addColor( DGtal::Color( 200, 10, 10 ) );    // red

	// lambda to convert a DGtal::Color object into a SVG::Color object
	auto DGtalColor2SVGColor = [](const DGtal::Color & c) { return SVG::Color(c.red(), c.green(), c.blue()); };

	std::vector< std::shared_ptr<SVG::AnimatedElement> > lines_ptr(segments.size(), nullptr);	// to each segment its line element

	// Starting from the end because the data is ordered as the CCO algo output them;
	// the last segment is guaranteed to be terminal, and the second to last is its brother (same parent segment)
	// rit is a reverse iterator to the last segment, which is created in the animation at this loop
	for(auto rit = segments.rbegin(); rit != segments.rend() && rit+1 != segments.rend(); rit += 2)
	{
		// find the reverse iterators and the indexes for the brother and the parent
		auto bro_rit = rit + 1;				// brother reverse iterator

		// parent reverse iterator
		unsigned int proxital_i = rit->myProxitalIndex;
		auto is_parent = [&proxital_i] (const Segment & s) { return s.myDistalIndex == proxital_i; };
		auto parent_rit = std::find_if(rit, segments.rend(), is_parent);

		// indexes 
		std::size_t i = segments.size() - 1 - (rit - segments.rbegin());
		std::size_t bro_i = segments.size() - 1 - (bro_rit - segments.rbegin());
		std::size_t parent_i = segments.size() - 1 - (parent_rit - segments.rbegin());		 

		// handy points for the next computations
		// the junction is the point where the brother, the parent and the current segment meet
		// the intersection is an arbitrary point for the sake of the animation; I think of it as the junction before the segment is born
		PointD<2> junction = myTree.myPoints[rit->myProxitalIndex];
		PointD<2> intersection;
		
		bool res = GeomHelpers::lineIntersection(myTree.myPoints[parent_rit->myProxitalIndex], myTree.myPoints[bro_rit->myDistalIndex],
									myTree.myPoints[rit->myProxitalIndex], myTree.myPoints[rit->myDistalIndex],
									intersection);
		if(!res)	// the lines are almost parallel or coincident (very unelikely)
		{
			// set intersection at starting point of added segment
			intersection = junction;
		}

		// define timestamps for this animation
		int start = duration - ((rit + 2 - segments.rbegin()) / 2) * segment_growth_dur;
		int step = start + delay;
		int end = start + segment_growth_dur;

		// initialize the lines with their color and thickness if they don't exist yet
		initializeSVGLine(myTree.myPoints[rit->myProxitalIndex],myTree.myPoints[rit->myDistalIndex], 
		   myTree.myRadii[rit->myDistalIndex], DGtalColor2SVGColor(cmap_grad(std::log(rit->myFlow))), duration,
		   lines_ptr[i]);

		initializeSVGLine(myTree.myPoints[bro_rit->myProxitalIndex],myTree.myPoints[bro_rit->myDistalIndex], 
		   myTree.myRadii[bro_rit->myDistalIndex], DGtalColor2SVGColor(cmap_grad(std::log(bro_rit->myFlow))), duration,
		   lines_ptr[bro_i]);

		initializeSVGLine(myTree.myPoints[parent_rit->myProxitalIndex],myTree.myPoints[parent_rit->myDistalIndex], 
		   myTree.myRadii[parent_rit->myDistalIndex], DGtalColor2SVGColor(cmap_grad(std::log(parent_rit->myFlow))), duration,
		   lines_ptr[parent_i]);

		// create the animations
		// added segment
		SVG::Animation * a_x1 = lines_ptr[i]->getAnimation("x1");
		a_x1->getTimelineRef().addKeyTime(step, intersection[0]);
		a_x1->getTimelineRef().addKeyTime(end, junction[0]);
		
		SVG::Animation * a_y1 = lines_ptr[i]->getAnimation("y1");
		a_y1->getTimelineRef().addKeyTime(step, intersection[1]);
		a_y1->getTimelineRef().addKeyTime(end, junction[1]);

		SVG::Animation * a_x2 = lines_ptr[i]->getAnimation("x2");
		a_x2->getTimelineRef().addKeyTime(step, intersection[0]);
		a_x2->getTimelineRef().addKeyTime(end, myTree.myPoints[rit->myDistalIndex][0]);

		SVG::Animation * a_y2 = lines_ptr[i]->getAnimation("y2");
		a_y2->getTimelineRef().addKeyTime(step, intersection[1]);
		a_y2->getTimelineRef().addKeyTime(end, myTree.myPoints[rit->myDistalIndex][1]);

		SVG::Animation * a_opacity = lines_ptr[i]->getAnimation("opacity");
		a_opacity->getTimelineRef().addKeyTime(start, 0.0);
		a_opacity->getTimelineRef().addKeyTime(step, 1.0);

		// brother segment
		SVG::Animation * a_bx1 = lines_ptr[bro_i]->getAnimation("x1");
		a_bx1->getTimelineRef().addKeyTime(start, myTree.myPoints[bro_rit->myDistalIndex][0]);
		a_bx1->getTimelineRef().addKeyTime(step, intersection[0]);
		a_bx1->getTimelineRef().addKeyTime(end, junction[0]);
		
		SVG::Animation * a_by1 = lines_ptr[bro_i]->getAnimation("y1");
		a_by1->getTimelineRef().addKeyTime(start, myTree.myPoints[bro_rit->myDistalIndex][1]);
		a_by1->getTimelineRef().addKeyTime(step, intersection[1]);
		a_by1->getTimelineRef().addKeyTime(end, junction[1]);

		SVG::Animation * a_bx2 = lines_ptr[bro_i]->getAnimation("x2");
		a_bx2->getTimelineRef().addKeyTime(start, myTree.myPoints[bro_rit->myDistalIndex][0]);

		SVG::Animation * a_by2 = lines_ptr[bro_i]->getAnimation("y2");
		a_by2->getTimelineRef().addKeyTime(start, myTree.myPoints[bro_rit->myDistalIndex][1]);

		SVG::Animation * a_bopacity = lines_ptr[bro_i]->getAnimation("opacity");
		a_bopacity->getTimelineRef().addKeyTime(start, 0.0);
		a_bopacity->getTimelineRef().addKeyTime(step, 1.0);

		// parent segment
		SVG::Animation * a_px2 = lines_ptr[parent_i]->getAnimation("x2");
		a_px2->getTimelineRef().addKeyTime(start, myTree.myPoints[bro_rit->myDistalIndex][0]);
		a_px2->getTimelineRef().addKeyTime(step, intersection[0]);
		a_px2->getTimelineRef().addKeyTime(end, junction[0]);

		SVG::Animation * a_py2 = lines_ptr[parent_i]->getAnimation("y2");
		a_py2->getTimelineRef().addKeyTime(start, myTree.myPoints[bro_rit->myDistalIndex][1]);
		a_py2->getTimelineRef().addKeyTime(step, intersection[1]);
		a_py2->getTimelineRef().addKeyTime(end, junction[1]);

		// edit parent distal
		parent_rit->myDistalIndex = bro_rit->myDistalIndex;
	}

	// adjust opacity of root
	SVG::Animation * a_opacity_root = lines_ptr[0]->getAnimation("opacity");
	a_opacity_root->getTimelineRef().addKeyTime(0.0, 1.0);

	// compute svg viewbox
	PointD<2> ub, lb;						// upper and lower bounds
	compBB<2>(myTree.myPoints, ub, lb);

	SVG::Svg svg(lb[0], lb[1],		// top left coordinates
			   ub[0] - lb[0],		// width
			   ub[1] - lb[1]);		// height

	// add elements to svg
	for(std::shared_ptr<SVG::AnimatedElement> & line_ptr : lines_ptr)
	{
		svg.addElement(std::move(line_ptr));
	}

	// write svg to file
	std::ofstream file;
	file.open(filename + ".svg");

	if(file.is_open())
	{
		file << svg;				// Write SVG to file

		file.close();	
	}
	else
	{
		std::cout << "Couldn't open " << filename  + ".svg" << std::endl;
	}
}



template<int TDim>
Image<TDim> TreeImageRenderer<TDim>::skeletonRender(unsigned int width)
{
	// initialize a Image<TDim> with the desired width and margin of 5%
	Image<TDim> organ_img(createDomainImage(width, width/20));
	Image<TDim> skeleton_render(organ_img.domain());		// image with the same dimensions

	// loop over segments, draw them with drawBresenhamLine
	for(const Segment & s : myTree.mySegments)
	{
		drawBresenhamLine<TDim>(skeleton_render,
			PointI<TDim>(myTree.myPoints[s.myProxitalIndex]),
			PointI<TDim>(myTree.myPoints[s.myDistalIndex]),
			255.0);
	}

	if(myOrganDomain.isDefined)
	{
		// blend the organ and the flow
		auto customblend = [](double a, double b) { return (b <= 0.001 ? a : b); };
		imageBlend<TDim>(organ_img, skeleton_render, customblend, skeleton_render);
	}

	return skeleton_render;
}



template<int TDim>
Image<TDim> TreeImageRenderer<TDim>::realisticRender(double sigma, unsigned int width)
{
	// initialize a Image<TDim> with a render of the flow of the artery tree
	Image<TDim> realistic_render(flowRender(width));

	// random
	std::random_device rd;
	std::mt19937 generator {rd()};

	// standard deviation multiplier when A == 0
	double c = std::sqrt(2 - M_PI/2);

	for(auto it = realistic_render.begin(); it != realistic_render.end(); it++)
	{
		// rician distribution parameters
		double A = *it;
		double s = (fabs(A) <= 0.001 ? sigma * c : sigma);

		std::normal_distribution<> nd_x(A, s);
		std::normal_distribution<> nd_y(0, s);

		*it = PointD<TDim>(nd_x(generator), nd_y(generator)).norm();
	}

	return realistic_render;
}



template<int TDim>
Image<TDim> TreeImageRenderer<TDim>::createDomainImage(unsigned int width, unsigned int margin_thickness)
{
	// if the organ domain is defined, no need to rescale the points
	if(myOrganDomain.isDefined)
	{

		Image<TDim> organ(myOrganDomain.myDomainMask);

		normalizeImageValues<TDim>(organ, 0.0, 128.0);		// organ domain has value 128.0

		return organ;
	}

	// Compute the coordinates of the bounding box containing all the points
	PointD<TDim> tree_lowerbound;
	PointD<TDim> tree_upperbound;

	compBB<TDim>(myTree.myPoints, tree_upperbound, tree_lowerbound);

	PointD<TDim> tree_size = tree_upperbound - tree_lowerbound;

	// factor between tree size and image size (without margins)
	double k = ((double) (width- 2*margin_thickness)) / tree_size[0];
	
	PointI<TDim> image_size;
	image_size[0] = width;
	for(unsigned int i = 1; i < TDim; i++)		// start at 1 because 0 is already done
	{
		image_size[i] = tree_size[i] * k + 2 * margin_thickness; 
	}

	// scale and move the points of the tree so their position reflect their position in the image
	// their position are still PointD<TDim> (real points)
	PointD<TDim> offset = PointD<TDim>::diagonal(margin_thickness) - k*tree_lowerbound; // + myDomain.lowerBound() but it's (0, 0)

	for(PointD<TDim> &p : myTree.myPoints)
	{
		p *= k;
		p += offset;
	}

	// scale the radii
	for(double &r : myTree.myRadii)
	{
		r *= k;
	}

	// return the Image<TDim> object
	// offset by one so that the size is valid)
	Image<TDim> res( Domain<TSpace>(PointI<TDim>(), image_size - PointI<TDim>::diagonal(1)) );
	normalizeImageValues<TDim>(res, 0.0, 0.0); 		// should be a black image after this

	return res; 
}



///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////               Other functions                 ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


template<>
void saveRender<2>(const Image<2> & image,
				const std::string & filename)
{
	auto min_val = std::min_element(image.constRange().begin(), image.constRange().end());
	auto max_val = std::max_element(image.constRange().begin(), image.constRange().end());

	if(*min_val == *max_val)
	{		// important because DGtal::GradientColorMap causes a segfault when min = max
		std::cout << "The image has only one value, no point displaying it.";
		return;
	}

	DGtal::GradientColorMap<double> gradient_cmap(*min_val, *max_val);

	gradient_cmap.addColor(DGtal::Color::Black);
	gradient_cmap.addColor(DGtal::Color::White);

	DGtal::STBWriter< Image<2>, DGtal::GradientColorMap<double> > 
		::exportPNG(filename + ".png", image, gradient_cmap);

	std::cout << "Render exported to " << filename << ".png" << std::endl;
}



template<>
void saveRender<3>(const Image<3> & image,
				const std::string & filename)
{
	Image<3> norm_image(image);
	normalizeImageValues<3>(norm_image, 0.0, 255.0);

	DGtal::functors::Cast<unsigned char> cast_functor;

	DGtal::VolWriter< Image<3>, DGtal::functors::Cast<unsigned char> >
		::exportVol(filename + ".vol", norm_image, true, cast_functor);

	// test nifti
	DGtal::ITKWriter< Image<3> >::exportITK(filename + ".nii", norm_image);


	std::cout << "Render exported to " << filename << ".nii and " << filename << ".vol" << std::endl;
}



template<int TDim>
void compBB(const std::vector< PointD<TDim> > & points,
			PointD<TDim> & upperbound, 
			PointD<TDim> & lowerbound)
{
	upperbound = points[0];
	lowerbound = points[0];

	for(auto &p : points)
	{
		for(unsigned int i = 0; i < TDim; i++)
		{
			if(p[i] < lowerbound[i])
			{
				lowerbound[i] = p[i];
			}
			else if(p[i] > upperbound[i])
			{
				upperbound[i] = p[i];
			}
		}
	}
}



template<int TDim>
bool projectOnStraightLine(const PointD<TDim>& ptA,
						   const PointD<TDim> & ptB,
						   const PointD<TDim> & ptC,
						   PointD<TDim> & ptP)
{
	if(ptA == ptB)
	{
		// there is no line if A == B
		throw std::runtime_error("projectOnStraightLine error : A and B are equals");
	}

    if (ptA == ptC || ptB == ptC)
    {
        ptP = ptC;
        return true;
    }

    PointD<TDim> vAB = ptB - ptA;
    PointD<TDim> vABn = vAB / vAB.norm();		// norm can't be 0

    PointD<TDim> vAC = ptC-ptA;
    double distPtA_Proj = vAC.dot(vABn);

    ptP = ptA + vABn * distPtA_Proj;

    PointD<TDim> vPA = ptA - ptP;
    PointD<TDim> vPB = ptB - ptP;
    
    return vPB.dot(vPA) <= 0 ;
}



bool initializeSVGLine(const PointD<2> & proxital,
					   const PointD<2> & distal, 
					   double radius,
					   const SVG::Color & color,
					   int duration,
					   std::shared_ptr<SVG::AnimatedElement> & line_ptr)
{
	if(!line_ptr)	// nullptr check
	{
		// make a shared pointer
		line_ptr = std::make_shared<SVG::Line>(proxital[0], proxital[1], distal[0], distal[1],radius * 3, color);

		// init the 5 attributes animation
		line_ptr->initializeAnimation("x1", proxital[0], proxital[0], duration, true, 1);
		line_ptr->initializeAnimation("y1", proxital[1], proxital[1], duration, true, 1);
		line_ptr->initializeAnimation("x2", distal[0], distal[0], duration, true, 1);
		line_ptr->initializeAnimation("y2", distal[1], distal[1], duration, true, 1);
		line_ptr->initializeAnimation("opacity", 0.0, 1.0, duration, false, 1);

		return true;
	}

	return false;
}


// explicit instantion (the application, at the moment, only supports 2D and 3D)
template class TreeImageRenderer<2>;
template class TreeImageRenderer<3>;