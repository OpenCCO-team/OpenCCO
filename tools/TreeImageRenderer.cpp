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
	/*
	auto segcmp = [&] (const Segment & s1, const Segment & s2)
	{
		return std::min(s1.myDistalIndex, s1.myProxitalIndex) < std::min(s2.myDistalIndex, s2.myProxitalIndex);
	};*/

	std::sort(mySegments.begin(), mySegments.end());
}



///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////            TreeImageRenderer                  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// Constructor

TreeImageRenderer::TreeImageRenderer(unsigned int width) 
	: myBackground( TDomain() ),
	myTreeImage( TDomain() ),
	myDistanceMap( TDomain() )
{
	importTreeData();
	myTree.sort();

	setImageSize(width, width/20);		// initializes myDomain
}



// Methods

void TreeImageRenderer::importTreeData()
{
	std::ifstream file;

	// get the vertex from vertex.dat
	file.open("vertex.dat");

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
		std::cout << "Couldn't open vertex.dat" << std::endl;
	}

	// get the edge data from edges.dat
	file.open("edges.dat");

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
		std::cout << "Couldn't open edges.dat" << std::endl;
	}

	// get the radii from from radius.dat
	file.open("radius.dat");

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
		std::cout << "Couldn't open radius.dat" << std::endl;
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



void TreeImageRenderer::createDistanceAlphaChannel()
{
	for(const TPoint &p : myDomain)
	{
		/*
		auto it = myTree.mySegments.begin();
		while (it != myTree.mySegments.end() && it->)
		*/
		return;
	}
}



// Other useful fonctions

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