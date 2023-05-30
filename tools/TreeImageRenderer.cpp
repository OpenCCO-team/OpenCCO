#include "TreeImageRenderer.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                   Segment                     ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// Operator overloading
bool operator<(const Segment & s1, const Segment & s2)
{
	return std::min(s1.myDistalIndex, s1.myProxitalIndex) < std::min(s2.myDistalIndex, s2.myProxitalIndex);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                 ArteryTree                    ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// Method

void ArteryTree::sort()
{
	// create the vector of the points' first coordinate offset by the point radii
	std::vector<double> p_offset(myRadii);
	for(unsigned int i = 0; i < p_offset.size(); i++)
	{
		p_offset[i] -= myPoints[i][0];
		p_offset[i] *= -1;
	}

	// sort a vector of indices (permutation)
	std::vector<std::size_t> p = sort_permutation(p_offset,
				[](const double & a, const double & b) { return a < b; });

	// apply permutation to myRadii, myPoints, and to the distal and proxital indices of the segments
	myRadii = apply_permutation(myRadii, p);
	myPoints = apply_permutation(myPoints, p);

	for(Segment & s : mySegments)
	{
		s.myDistalIndex = p[s.myDistalIndex];
		s.myProxitalIndex = p[s.myProxitalIndex];
	}

	// finally sort the mySegment vector, based on the minimum of myProxitalIndex and myDistalIndex
	std::sort(mySegments.begin(), mySegments.end());
}



///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////            TreeImageRenderer                  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// Constructor

TreeImageRenderer::TreeImageRenderer(const unsigned int width,
									 const std::string & radii_filename,
									 const std::string & vertices_filename,
									 const std::string & edges_filename)
	: myBackground( TDomain() ),
	myTreeImage( TDomain() ),
	myDistanceMap( TDomain() )
{
	importTreeData(radii_filename, vertices_filename, edges_filename);
	//myTree.sort();

	setImageSize(width, width/20);		// initializes myDomain
}



// Methods

void TreeImageRenderer::importTreeData(const std::string & radii_filename,
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
			iss >> p[0] >> p[1];
			
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



void TreeImageRenderer::setImageSize(unsigned int width, unsigned int margin_thickness)
{
	// Compute the coordinates of the bounding box containing all the points
	TPointD tree_lowerbound(myTree.myPoints[0]);
	TPointD tree_upperbound(myTree.myPoints[0]);

	compBB(myTree.myPoints, tree_upperbound, tree_lowerbound);

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
	myBackground = TImage2D(myDomain);
	myTreeImage = TImage2D(myDomain);
	myDistanceMap = TImage2D(myDomain);

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

	compBB(myTree.myPoints, tree_upperbound, tree_lowerbound);

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



void TreeImageRenderer::createDistanceMap()
{
	for(const TPoint &p : myDomain)
	{
		// base value is 0, since we're searching for minimal distances it would stop the algorithm instantly
		// infinity ensure that any distance we check will be accepted
		myDistanceMap.setValue(p, std::numeric_limits<double>::infinity());

		auto it = myTree.mySegments.begin();

		// the while loop will iterate over the segment in their order of appearance from left to right
		// like a sweeping line (segment must be sorted, it's taken care of in the constructor of TreeImageRenderer)
		std::size_t ind_point = std::min(it->myDistalIndex, it->myProxitalIndex);	// leftmost point
		double sweeping_line_x = myTree.myPoints[ind_point][0] - myTree.myRadii[ind_point];

		while (it != myTree.mySegments.end()
			/*&& myDistanceMap(p) > sweeping_line_x - p[0]*/)		// check if the sweeping line is too far compared to the already register min distance 
		{
			// project p onto the segment
			TPointD proj;
			bool isproj = projectOnStraightLine(myTree.myPoints[it->myDistalIndex],
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




void TreeImageRenderer::createTreeImage()
{
	for(const TPoint &p : myDomain)
	{
		auto it = myTree.mySegments.begin();
		while (it != myTree.mySegments.end())
		{
			// project p onto the segment
			TPointD proj;
			bool isproj = projectOnStraightLine(myTree.myPoints[it->myDistalIndex],
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
					myTreeImage.setValue(p, 1.0);
				}
			}

			// update 
			it++;
		}
	}
}

const TImage2D & TreeImageRenderer::distanceMap() const
{
	return myDistanceMap;
}

const TImage2D & TreeImageRenderer::treeImage() const
{
	return myTreeImage;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////               Other functions                 ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////



void compBB(std::vector<TPointD> &points, TPointD &upperbound, TPointD &lowerbound)
{
	for(TPointD &p : points)
	{
		for(unsigned int i = 0; i < TPointD::dimension; i++)
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




bool projectOnStraightLine(const TPointD & ptA,
						   const TPointD & ptB,
						   const TPointD & ptC,
						   TPointD & ptP)
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

    TPointD vAB = ptB - ptA;
    TPointD vABn = vAB / vAB.norm();		// norm can't be 0

    TPointD vAC = ptC-ptA;
    double distPtA_Proj = vAC.dot(vABn);

    ptP = ptA + vABn * distPtA_Proj;

    TPointD vPA = ptA - ptP;
    TPointD vPB = ptB - ptP;
    
    return vPB.dot(vPA) <= 0 ;
}



TDomain domainIntersect(TDomain &dom1, TDomain &dom2)
{
	TPoint upperbound;
	TPoint lowerbound;

	for(unsigned int i = 0; i < TPoint::dimension; ++i)
	{
		lowerbound[i] = std::max(dom1.lowerBound()[i], dom2.lowerBound()[i]);
		upperbound[i] = std::min(dom1.upperBound()[i], dom2.upperBound()[i]);
	}

	return TDomain(lowerbound, upperbound);
}