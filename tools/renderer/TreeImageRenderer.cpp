#include "TreeImageRenderer.h"


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/io/colormaps/GradientColorMap.h"

#include "DGtal/io/writers/PGMWriter.h"
#include "DGtal/io/writers/STBWriter.h"
#include "DGtal/io/writers/VolWriter.h"


///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////          TreeImageRenderer<TDim>              ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// Constructor

template<int TDim>
TreeImageRenderer<TDim>::TreeImageRenderer(const unsigned int width,
									 const std::string & radii_filename,
									 const std::string & vertices_filename,
									 const std::string & edges_filename)
	: myBackground( TDomain() ),
	myTreeImage( TDomain() ),
	myDistanceMap( TDomain() )
{
	importTreeData(radii_filename, vertices_filename, edges_filename);

	setImageSize(width, width/20);		// initializes myDomain
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
	/*
	DGtal::VolWriter< TImage, DGtal::functors::Cast<unsigned char> > 
			::exportVol(filename + ".vol", myTreeImage, true, cast_functor);
	*/

	DGtal::VolWriter<TImage, DGtal::functors::Cast<unsigned char> >
		::exportVol(filename + ".vol", myTreeImage, true, cast_functor);
}



///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////               Other functions                 ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////



template<int TDim>
void compBB(std::vector< typename TreeImageRenderer<TDim>::TPointD > &points,
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

// explicit instantion (the application, at the moment, only supports 2D and 3D)
template class TreeImageRenderer<2>;
template class TreeImageRenderer<3>;