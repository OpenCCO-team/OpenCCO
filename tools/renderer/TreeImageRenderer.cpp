#include "TreeImageRenderer.h"

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/PGMWriter.h"
#include "DGtal/io/writers/STBWriter.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/io/writers/ITKWriter.h"

#include "GeomHelpers.h"

#include "svg_elements.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////          TreeImageRenderer<TDim>              ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// Constructor

template<int TDim>
TreeImageRenderer<TDim>::TreeImageRenderer(const std::string & radii_filename,
									 const std::string & vertices_filename,
									 const std::string & edges_filename)
	: myBackground( TDomain() ), myTreeImage( TDomain() ), myDistanceMap( TDomain() )
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
			TPointD p;

			std::istringstream iss(line);

			for(unsigned int i = 0; i < TPointD::dimension; i++)
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

	std::cout << "Imported " << myTree.myPoints.size() << " vertices." << std::endl;
	std::cout << "Imported " << myTree.mySegments.size() << " edges." << std::endl;
	std::cout << "Imported " << myTree.myRadii.size() << " radii." << std::endl;
}



template<int TDim>
void TreeImageRenderer<TDim>::setImageSize(unsigned int width, unsigned int margin_thickness)
{
	// Compute the coordinates of the bounding box containing all the points
	TPointD tree_lowerbound(myTree.myPoints[0]);
	TPointD tree_upperbound(myTree.myPoints[0]);

	compBB<TDim>(myTree.myPoints, tree_upperbound, tree_lowerbound);

	TPointD tree_size = tree_upperbound - tree_lowerbound;

	// factor between tree size and image size (without margins)
	double k = ((double) (width- 2*margin_thickness)) / tree_size[0];
	
	TPoint image_size;
	image_size[0] = width;
	for(unsigned int i = 1; i < TPoint::dimension; i++)		// start at 1 because 0 is already done
	{
		image_size[i] = tree_size[i] * k + 2 * margin_thickness; 
	}

	myDomain = TDomain(TPoint(), image_size - TPoint::diagonal(1)); // offset by one so that the size is valid

	// set the size of the images
	myBackground = TImage(myDomain);
	myTreeImage = TImage(myDomain);
	myDistanceMap = TImage(myDomain);

	// scale and move the points of the tree so their position reflect their position in the image
	// their position are still TPointD (real points)
	TPointD offset = TPointD::diagonal(margin_thickness) - k*tree_lowerbound; // + myDomain.lowerBound() but it's (0, 0)

	for(TPointD &p : myTree.myPoints)
	{
		p *= k;
		p += offset;
	}

	// scale the radii
	for(double &r : myTree.myRadii)
	{
		r *= k;
	}

	// verbose debugging :)

	std::cout << "Image size : ";
	for(int &a : image_size)
	{
		std::cout << a << " ";
	}
	std::cout << std::endl;

	tree_lowerbound = myTree.myPoints[0];
	tree_upperbound = myTree.myPoints[0];

	compBB<TDim>(myTree.myPoints, tree_upperbound, tree_lowerbound);

	std::cout << "Tree bounds : " << std::endl;
	for (unsigned int i = 0; i < TPointD::dimension; i++)
	{
		std::cout << tree_lowerbound[i] << " ";
	}
	std::cout << std::endl;
	for (unsigned int i = 0; i < TPointD::dimension; i++)
	{
		std::cout << tree_upperbound[i] << " ";
	}
	std::cout << std::endl;
}



template<int TDim>
void TreeImageRenderer<TDim>::createDistanceMap()
{
	for(const TPoint &p : myDomain)
	{
		// base value is 0, since we're searching for minimal distances it would stop the algorithm instantly
		// infinity ensure that any distance we check will be accepted
		myDistanceMap.setValue(p, std::numeric_limits<double>::infinity());

		auto it = myTree.mySegments.begin();

		// the while loop will iterate over the segment in their order of appearance from left to right
		// like a sweeping line (segment must be sorted, it's taken care of in the constructor of TreeImageRenderer<TDim>)
		std::size_t ind_point = std::min(it->myDistalIndex, it->myProxitalIndex);	// leftmost point
		double sweeping_line_x = myTree.myPoints[ind_point][0] - myTree.myRadii[ind_point];

		while (it != myTree.mySegments.end()
			/*&& myDistanceMap(p) > sweeping_line_x - p[0]*/)		// check if the sweeping line is too far compared to the already register min distance 
		{
			// project p onto the segment
			TPointD proj;
			bool isproj = projectOnStraightLine<TDim>(myTree.myPoints[it->myDistalIndex],
												myTree.myPoints[it->myProxitalIndex],
												p,
												proj);

			if(isproj)		// the projection belongs to the segment
			{
				double dist_pproj = (proj - p).norm();

				double dist_proxitalproj = (proj - myTree.myPoints[it->myProxitalIndex]).norm();
				double dist_proxitaldistal = (myTree.myPoints[it->myDistalIndex] - myTree.myPoints[it->myProxitalIndex]).norm();
				double interpolated_radius = (myTree.myRadii[it->myDistalIndex] - myTree.myRadii[it->myProxitalIndex]) * dist_proxitalproj/dist_proxitaldistal
											+ myTree.myRadii[it->myProxitalIndex];

				if(myDistanceMap(p) > dist_pproj - interpolated_radius)	
				{
					myDistanceMap.setValue(p, std::max(dist_pproj - interpolated_radius, 0.0));
				}
			}
			else		// the projection doesn't belong to the segment
			{			// the distance from p to the segment is either the distance to the distal or the distance to the proxital
						// offset by their respective radius
				double dist_pproxital = (myTree.myPoints[it->myProxitalIndex] - p).norm() - myTree.myRadii[it->myProxitalIndex];
				double dist_pdistal = (myTree.myPoints[it->myDistalIndex] - p).norm() - myTree.myRadii[it->myDistalIndex];

				double min_dist = std::min(dist_pproxital, dist_pdistal);

				if(myDistanceMap(p) > min_dist)
				{
					myDistanceMap.setValue(p, min_dist);
				}
			}

			// update 
			it++;
			if( it != myTree.mySegments.end() )
			{
				ind_point = std::min(it->myDistalIndex, it->myProxitalIndex);
				sweeping_line_x = myTree.myPoints[ind_point][0] - myTree.myRadii[ind_point];	
			}		
		}
	}
}




template<int TDim>
void TreeImageRenderer<TDim>::createTreeImage()
{
	for(const TPoint &p : myDomain)
	{
		myTreeImage.setValue(p, 255);

		auto it = myTree.mySegments.begin();
		while (it != myTree.mySegments.end())
		{
			// project p onto the segment
			TPointD proj;
			bool isproj = projectOnStraightLine<TDim>(myTree.myPoints[it->myDistalIndex],
												myTree.myPoints[it->myProxitalIndex],
												p,
												proj);

			if(isproj)		// the projection belongs to the segment
			{
				double dist_pproj = (proj - p).norm();

				double dist_proxitalproj = (proj - myTree.myPoints[it->myProxitalIndex]).norm();
				double dist_proxitaldistal = (myTree.myPoints[it->myDistalIndex] - myTree.myPoints[it->myProxitalIndex]).norm();
				double interpolated_radius = (myTree.myRadii[it->myDistalIndex] - myTree.myRadii[it->myProxitalIndex]) * dist_proxitalproj/dist_proxitaldistal
											+ myTree.myRadii[it->myProxitalIndex];

				if(dist_pproj < interpolated_radius)	
				{
					myTreeImage.setValue(p, 1.0);
				}
			}
			else		// the projection doesn't belong to the segment
			{			// the distance from p to the segment is either the distance to the distal or the distance to the proxital
						// offset by their respective radius
				double dist_pproxital = (myTree.myPoints[it->myProxitalIndex] - p).norm() - myTree.myRadii[it->myProxitalIndex];
				double dist_pdistal = (myTree.myPoints[it->myDistalIndex] - p).norm() - myTree.myRadii[it->myDistalIndex];

				double min_dist = std::min(dist_pproxital, dist_pdistal);

				if(min_dist < 0)
				{
					myTreeImage.setValue(p, 0);
				}
			}

			// update 
			it++;
		}
	}
}


template<>
void TreeImageRenderer<2>::saveRender(const std::string & filename)
{
	auto min_val = std::min_element(myTreeImage.constRange().begin(), myTreeImage.constRange().end());
	auto max_val = std::max_element(myTreeImage.constRange().begin(), myTreeImage.constRange().end());

	DGtal::GradientColorMap<double> gradient_cmap(*min_val, *max_val);

	gradient_cmap.addColor(DGtal::Color::Black);
	gradient_cmap.addColor(DGtal::Color::White);

	DGtal::STBWriter< TImage, DGtal::GradientColorMap<double> > 
		::exportPNG(filename + ".png", myTreeImage, gradient_cmap);
}



template<>
void TreeImageRenderer<3>::saveRender(const std::string & filename)
{
	
	DGtal::functors::Cast<unsigned char> cast_functor;

	DGtal::VolWriter< TImage, DGtal::functors::Cast<unsigned char> >
		::exportVol(filename + ".vol", myTreeImage, true, cast_functor);

	// test nifti
	DGtal::ITKWriter<TImage>::exportITK("render.nii", myTreeImage);
}


template<>
void TreeImageRenderer<2>::treeConstructionAnimation(const std::string & filename, int duration)
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

		// handy points for what follows
		// the junction is the point where the brother, the parent and the current segment meet
		// the intersection is an arbitrary point for the sake of the animation; I think of it as the junction before the segment is born
		TPointD junction = myTree.myPoints[rit->myProxitalIndex];
		TPointD intersection;
		
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
	TPointD ub, lb;						// upper and lower bounds
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
	file.open(filename);

	if(file.is_open())
	{
		file << svg;				// Write SVG to file

		file.close();	
	}
	else
	{
		std::cout << "Couldn't open " << filename << std::endl;
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////               Other functions                 ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////



template<int TDim>
void compBB(const std::vector< typename TreeImageRenderer<TDim>::TPointD > &points,
			typename TreeImageRenderer<TDim>::TPointD &upperbound, 
			typename TreeImageRenderer<TDim>::TPointD &lowerbound)
{
	for(auto &p : points)
	{
		for(unsigned int i = 0; i < TreeImageRenderer<TDim>::myDim; i++)
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
bool projectOnStraightLine(const typename TreeImageRenderer<TDim>::TPointD & ptA,
						   const typename TreeImageRenderer<TDim>::TPointD & ptB,
						   const typename TreeImageRenderer<TDim>::TPointD & ptC,
						   typename TreeImageRenderer<TDim>::TPointD & ptP)
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

    auto vAB = ptB - ptA;
    auto vABn = vAB / vAB.norm();		// norm can't be 0

    auto vAC = ptC-ptA;
    double distPtA_Proj = vAC.dot(vABn);

    ptP = ptA + vABn * distPtA_Proj;

    auto vPA = ptA - ptP;
    auto vPB = ptB - ptP;
    
    return vPB.dot(vPA) <= 0 ;
}



bool initializeSVGLine(const TreeImageRenderer<2>::TPointD & proxital,
					   const TreeImageRenderer<2>::TPointD & distal, 
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