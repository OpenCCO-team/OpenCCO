/**
 *  OpenCCO implementation 
 *  Copyright (C) 2023 B. Kerautret;  Phuc Ngo, N. Passat H. Talbot and C. Jaquet
 *
 *  This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 **/


#pragma once

#if defined(GEOMHELPERS_RECURSES)
#error Recursive header files inclusion detected in Geomhelpers.h
#else // defined(GEOMHELPERS_RECURSES)
/** Prevents recursive inclusion of headers. */
#define GEOMHELPERS_RECURSES

#if !defined GEOMHELPERS_h
/** Prevents repeated inclusion of headers. */
#define GEOMHELPERS_h
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/shapes/implicit/ImplicitBall.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/shapes/EuclideanShapesDecorator.h"



#include "ceres/ceres.h"
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;

// types used through the algorithm

template <int TDim>
using PointI = DGtal::PointVector<TDim, int>;

template <int TDim>
using PointD = DGtal::PointVector<TDim, double>;

template <class TPoint>
using Space = DGtal::SpaceND<TPoint::dimension, typename TPoint::Component>;

template <class Space>
using Domain = DGtal::HyperRectDomain<Space>;


namespace GeomHelpers {




/**
 * Computes the projection of a Point \a ptC on the real line defined by the two points (\a ptA, \a ptB), and
 * return true if the projected point is inside the segment closed interval [A,B].
 *
 * @param[in] ptA one of the two points defining the straight line.
 * @param[in] ptB one of the two points defining the straight line.
 * @param[in] ptC the point to be projected.
 * @param[out] ptProjected the projected point.
 * @return true if ptProjected is inside the segment [A,B].
 **/

template<int TDim>
inline
bool
projectOnStraightLine(const PointD<TDim> & ptA,
					  const PointD<TDim> & ptB,
					  const PointD<TDim> & ptC,
					  PointD<TDim> & ptProjected)
{
	if (ptA==ptC)
	{
		for(auto i=0; i<TDim; i++){ ptProjected[i]=ptA[i];}
		return true;
	}
	if (ptB==ptC)
	{
		for(auto i=0; i<TDim; i++){ ptProjected[i]=ptB[i];}
		return true ;
	}

	PointD<TDim> vAB = ptB - ptA;//(ptB[0]- ptA[0], ptB[1]- ptA[1]);
	PointD<TDim> vABn = vAB / vAB.norm();//((double)vAB[0], (double)vAB[1]);
	//double norm = vABn.norm();
	//vABn[0] /= norm;
	//vABn[1] /= norm;

	PointD<TDim> vAC = ptC-ptA;// (ptC[0]-ptA[0], ptC[1]-ptA[1]);
	double distPtA_Proj = vAC.dot(vABn);

	for(std::size_t i=0; i<TDim; i++){ ptProjected[i]= ptA[i]+vABn[i]*(distPtA_Proj);}
	//ptProjected[0] = ptA[0]+vABn[0]*(distPtA_Proj);
	//ptProjected[1] = ptA[1]+vABn[1]*(distPtA_Proj);
	bool res = false;
	for(std::size_t i=0; i<TDim; i++) { res = res || (ptA[i]<ptB[i] && ptProjected[i]<=ptB[i]); }
	return distPtA_Proj>=0 && res;
	/*
	return  distPtA_Proj>=0 && ((ptA[0]<ptB[0] && ptProjected[0]<=ptB[0] ) ||
	(ptA[0]>ptB[0] && ptProjected[0]>=ptB[0] ) ||
	(ptA[0]==ptB[0] && ptA[1]<ptB[1] && ptProjected[1]<=ptB[1]) ||
	(ptA[0]==ptB[0] && ptA[1]>=ptB[1] && ptProjected[1]>=ptB[1]));
	*/
}




/**
 * From two segments represented by the end points, it returns their
 * Euclidean distance
 * @param segA first point of the first segment.
 * @param segB second point of the first segment.
 * @param segC first point of the first segment.
 * @param segD second point of the first segment.
 * @return distance between two segments
 **/
template<int TDim>
inline
double
segment2segmentDistance(const PointD<TDim> &segA, const PointD<TDim> &segB,
				const PointD<TDim> &segC, const PointD<TDim> &segD)
{
	PointD<TDim> A, B, C, D;
	PointD<TDim> d1, d2, d12;
	bool sameAB = true;
	for (int i = 0; i<TDim; i++) {
		A[i] = segA[i];
		B[i] = segB[i];
		C[i] = segC[i];
		D[i] = segD[i];

		d1[i] = segB[i] - segA[i]; //d2[i] : di2
		if(abs(d1[i])!=0)
			sameAB = false;
	}

	if(sameAB == true) { //switch AB <-> CD
		A = segC;
		B = segD;
		C = segA;
		D = segB;
	}

	double D1 = 0, D2 = 0, R = 0, S1 = 0, S2 = 0, den = 0;
	bool sAB = true, sCD = true;
	for (int i = 0; i<TDim; i++) {
		d1[i] = B[i] - A[i]; //d1[i] : di1
		d2[i] = D[i] - C[i]; //d2[i] : di2
		d12[i] = C[i] - A[i]; //d12[i] : di12

		D1 += d1[i]*d1[i];
		D2 += d2[i]*d2[i];
		R += d1[i]*d2[i];
		S1 += d1[i]*d12[i];
		S2 += d2[i]*d12[i];

		if(abs(d1[i])!=0)
			sAB = false;
		if(abs(d2[i])!=0)
			sCD = false;
	}

	den = D1*D2 - R*R;

	double u = 0, t = 0, DD = 0;
	int step = 0;
	//Case 1 : degenerate case of A=B but C!=D
	if(sAB == false && sCD == true) {
		u = 0;
		step = 4;
	}
	//Case 2 : degenerate case of A=B and C=D
	else if(sAB == true && sCD == true) {
		u = 0;
		t = 0;
		step = 5;
	}
	//Case 3 : degenerate case of A!=B and C!=D and denominator = 0
	else if(sAB == false && sCD == false && den == 0) {
		t = 0;
		step = 3;
	}
	//Case 4 : otherwise
	else {
		step = 2;
	}

	//Process of each case
	if(step <=2 ) {
		t = (S1*D2 - S2*R)/den;
		if(t<0)
		t = 0;
		else if(t>1)
		t = 1;
	}

	if(step <=3 ) {
		u = (t*R - S2)/D2;
		if(u<0)
			u = 0;
		else {
		if(u>1)
			u = 1;
		else
			step = 5;
		}
	}

	if(step <= 4) {
		t = (u*R + S1)/D1;
		if(t<0)
			t = 0;
		else if(t>1)
			t = 1;
	}
	if(step <= 5) {
		double d=0;
		for(int i=0; i<TDim; i++) {
			d = (d1[i]*t - d2[i]*u - d12[i]);
			DD += d*d;
		}
	}
	return sqrt(DD);
}

template<int TDim>
inline
double
point2segmentDistance(const PointD<TDim> &segA, const PointD<TDim> &segB,
					  const PointD<TDim> &ptP)
{
	return segment2segmentDistance<TDim>(segA, segB, ptP, ptP);
}

/**
 * From two thick segments represented by their end points and radius,
 * it indicate if they intersect
 * @param segA first point of the first segment.
 * @param segB second point of the first segment.
 * @param rAB radius of the first segment.
 * @param segC first point of the first segment.
 * @param segD second point of the first segment.
 * @param rCD radius of the second segment.
 * @return bool indicating if two segments intersect
 **/
template<int TDim>
inline
bool
isIntersecting(const PointD<TDim> &segA, const PointD<TDim> &segB, double rAB,
			   const PointD<TDim> &segC, const PointD<TDim> &segD, double rCD)
{
	PointD<TDim> cAB = (segA + segB) / 2.0;
	PointD<TDim> cCD = (segC + segD) / 2.0;
	double lAB = (segA - segB).norm();
	double lCD = (segC - segD).norm();
	double dC = (cAB - cCD).norm();
	bool res = dC > ((lAB + lCD) / 2.0 + rAB + rCD);

	if(res==true)
		return false;

	double distanceSeg = segment2segmentDistance<TDim>(segA, segB, segC, segD);

	if(distanceSeg > (rAB + rCD))
		return false;

	return true;
}


/**
 * Calculates the intersection of 2 lines each defined by 2 points.
 * Only for 2D lines.
 * @param segA first point of the first line.
 * @param segB second point of the first line.
 * @param segC first point of the first line.
 * @param segD second point of the first line.
 * @param[out] The point where the two lines intersect.
 * @return bool indicating if two line intersect.
 * @note the formula comes from https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Given_two_points_on_each_line
 **/
inline
bool
lineIntersection(const PointD<2> &segA,
				 const PointD<2> &segB,
				 const PointD<2> &segC,
				 const PointD<2> &segD,
				 PointD<2> &intersection)
{
	double denominator = (segA[0] - segB[0]) * (segC[1] - segD[1]) - (segA[1] - segB[1]) * (segC[0] - segD[0]);

	if(std::fabs(denominator) < 0.0001)  // denominator = 0 implies that the lines are parallel or coincident
	{
		return false;
	}

	// variables for the next computation
	double v1 = segA[0] * segB[1] - segA[1] * segB[0];
	double v2 = segC[0] * segD[1] - segC[1] * segD[0];

	intersection[0] = v1 * (segC[0] - segD[0]) - v2 * (segA[0] - segB[0]);
	intersection[1] = v1 * (segC[1] - segD[1]) - v2 * (segA[1] - segB[1]);

	intersection /= denominator;

	return true;
}



template<typename TPoint, typename TDSet>
inline
TDSet
pointsOnSphere(const TPoint & ptCenter,
			   double radius)
{

	typedef DGtal::ImplicitBall< Space<TPoint> > MyBall;
	typedef DGtal::EuclideanShapesCSG< MyBall, MyBall > Minus;
	typedef DGtal::GaussDigitizer< Space<TPoint>, Minus > MyGaussDigitizer;

	MyBall disk( ptCenter, radius-0.5 );
	MyBall diskDilate( ptCenter, radius+0.5 );
	Minus border ( diskDilate );
	border.minus( disk );
	MyGaussDigitizer digShape;
	digShape.attach( border );
	digShape.init( border.getLowerBound(), border.getUpperBound(), 1 );
	Domain< Space<TPoint> > domainShape = digShape.getDomain();
	TDSet  aSet(domainShape);
	DGtal::Shapes<Domain< Space<TPoint> >>::digitalShaper( aSet, digShape );

	return aSet;
}

/**
 * From two segments represented by the end points, it returns true if
 * there exists an intersection and 0 in the other cases.
 *
 * @param seg1ptA first point of the first segment.
 * @param seg1ptB second point of the first segment.
 * @param seg2ptA first point of the first segment.
 * @param seg2ptB second point of the first segment.
 * @return true if an intersection point exist and false in other cases.
 * @note Note that the intersection point should be insides the segments.
 * If the lines supported by the segments are coincident or the line are parellel  the function returns true.
 **/


template <typename TPoint>
inline
bool
hasIntersection(const TPoint &seg1ptA, const TPoint &seg1ptB,
				const TPoint &seg2ptA, const TPoint &seg2ptB)
{
	double  d = ((seg2ptB[1] - seg2ptA[1])*(seg1ptB[0] - seg1ptA[0])) -
		((seg2ptB[0] - seg2ptA[0])*(seg1ptB[1] - seg1ptA[1]));
	double a = ((seg2ptB[0] - seg2ptA[0])*(seg1ptA[1] - seg2ptA[1])) -
		((seg2ptB[1] - seg2ptA[1])*(seg1ptA[0] - seg2ptA[0]));
	double b = ((seg1ptB[0] - seg1ptA[0])*(seg1ptA[1] - seg2ptA[1])) -
		((seg1ptB[1] - seg1ptA[1])*(seg1ptA[0] - seg2ptA[0]));

	if ( d==0.0 ) {
		// test coincident
		if (a==0.0 && b == 0.0 ) {
			return false;
		}
		else
		return false;
	}

	double ua = a / d;
	double ub = b / d;
	return ua > 0.0f && ua < 1.0f && ub > 0.0f && ub < 1.0f;
	// Get the intersection point.
	//intersection.x_ = begin_.x_ + ua*(end_.x_ - begin_.x_);
	//intersection.y_ = begin_.y_ + ua*(end_.y_ - begin_.y_);
}

struct CostOptPos {
	double deltap1;
	double deltap2;
	double f0;
	double f1;
	double f2;
	double l0;
	double l1;
	double l2;
	double gamma;

	template <typename T>
	bool operator()(const T* const x1, const T* const x2, T* residual) const {
		//residual[0] =  deltap1*x1[0]*x1[0]*pow((f0*(((x1[0]*x1[0]*x1[0])/f1)+(x2[0]*x2[0]*x2[0])/f2)), 2.0/3.0)
		//- (f0*l0*x1[0]*x1[0]) - f1 * l1 * pow(f0*((x1[0]*x1[0]*x1[0]/f1)+(x2[0]*x2[0]*x2[0])/f2), 2.0/3.0);
		//residual[1] =  deltap2*x2[0]*x2[0]*pow((f0*(((x1[0]*x1[0]*x1[0])/f1)+(x2[0]*x2[0]*x2[0])/f2)), 2.0/3.0)
		//- (f0*l0*x2[0]*x2[0]) - f2 * l2 * pow(f0*((x1[0]*x1[0]*x1[0]/f1)+(x2[0]*x2[0]*x2[0])/f2), 2.0/3.0);
		//Equation 28
		residual[0] =  deltap1*x1[0]*x1[0]*(pow(f0*((pow(x1[0],(3.0+gamma)/2.0)/f1) + (pow(x2[0],(3.0+gamma)/2.0)/f2)),2.0/3.0))
		- (f0*l0*x1[0]*x1[0]) - f1 * l1 * (pow(f0*((pow(x1[0],(3.0+gamma)/2.0)/f1) + (pow(x2[0],(3.0+gamma)/2.0))/f2),2.0/3.0));
		residual[1] =  deltap2*x2[0]*x2[0]*(pow(f0*((pow(x1[0],(3.0+gamma)/2.0)/f1) + (pow(x2[0],(3.0+gamma)/2.0)/f2)),2.0/3.0))
		- (f0*l0*x2[0]*x2[0]) - f2 * l2 * (pow(f0*((pow(x1[0],(3.0+gamma)/2.0)/f1) + (pow(x2[0],(3.0+gamma)/2.0))/f2),2.0/3.0));
		return true;
	}
};


/**
 * @return true if a solution exists
 */
static bool kamyiaOpt(double gamma, double deltaP1, double deltaP2, double f0, double f1, double f2, double l0, double l1, double l2, double &xx1, double &xx2) 
{
	CostOptPos *f = new CostOptPos();
	f->deltap1 = deltaP1;
	f->deltap2 = deltaP2;
	f->f0 = f0;
	f->f1 = f1;
	f->f2 = f2;
	f->l0 = l0;
	f->l1 = l1;
	f->l2 = l2;
	f->gamma = gamma;

	// const double initial_x = x;
	// Build the problem.
	Problem problem;
	// Set up the only cost function (also known as residual). This uses
	// auto-differentiation to obtain the derivative (jacobian).
	CostFunction* cost_function = new AutoDiffCostFunction<CostOptPos, 2, 1, 1>(f);
	problem.AddResidualBlock(cost_function, nullptr, &xx1, &xx2);
	// Run the solver!
	Solver::Options options;
	options.minimizer_progress_to_stdout = false;
	Solver::Summary summary;
	Solve(options, &problem, &summary);
	//std::cout << summary.BriefReport() << "\n";

	return summary.IsSolutionUsable();
}



template< typename TImage, typename TImageDistance>
inline
TImageDistance
getImageDistance(const TImage &image, unsigned int threshold=128)
{
	typedef typename TImage::Domain::Space ImageSpace;
	typedef DGtal::ExactPredicateLpSeparableMetric<ImageSpace, TImage::Domain::dimension> L2Metric;

	TImageDistance res (image.domain());
	typedef DGtal::functors::IntervalForegroundPredicate<TImage> Binarizer;
	typedef DGtal::DistanceTransformation<ImageSpace, Binarizer, L2Metric> DTL;
	L2Metric l2metric;
	Binarizer b(image, threshold, 255);
	DTL dt(&image.domain(),&b, &l2metric);
	for (auto p: dt.domain()){
		res.setValue(p, dt(p));
	}
	return res;
}


/**
 * @brief returns a list of points of a n-D sphere.
 * @brief The n-D sphere is of radius 1, centered at the origin.
 * @param n The amount of points to spread.
 * @returns a std::vector of n-D points.
 **/
template <int TDim>
std::vector< PointD<TDim> > evenlySpreadPoints(unsigned int n)
{
	std::vector< PointD<TDim> > points;
	if(TDim == 2)
	{
		const double delta_theta = 2 * M_PI / n;

		for(std::size_t i = 0; i < n; i++)
		{
			// angle
			double theta = i * delta_theta;
			points.emplace_back(cos(theta), sin(theta));
		}
	}
	else if(TDim == 3)
	{
		const double epsilon = 0.36;				// offset for the indices, proven best for most n values : https://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/#more-3069
		const double double_GR = (1 + sqrt(5));		// golden ratio times 2 (to save a multiplication by 2 later)

		for(std::size_t i = 0; i < n; i++)
		{
			// angles
			double phi = acos(1 - 2 * (i + epsilon) / n);
			double theta = M_PI * double_GR * i;

			points.emplace_back(cos(theta) * sin(phi),
							 sin(theta) * sin(phi),
							 cos(phi) );
		}
	}

	return points;
}

/**
 * @brief Computes a normalized orthogonal vector by solving the dot product being zero.
 * @brief In dimensions >= 2, there's more than one solution, this function arbitrarily chooses one with
 * @brief as much coordinates being zero as possible.
 * @param vec The vector for which we want an orthogonal vector.
 * @returns a PointD<2> vector with at least one non-zero coordinate (unless the input vector is full of zeros).
 **/
template<int TDim>
PointD<TDim> orthogonalVector(const PointD<TDim> & vec);

template <>
inline
PointD<2> orthogonalVector<2>(const PointD<2> & vec)
{
	double epsilon = 0.00001;
	PointD<2> vec_ortho(0.0, 0.0);

	if(fabs(vec[0]) < epsilon)	// too close to zero
	{
		if(fabs(vec[1]) < epsilon)
		{
			// (0.0, 0.0)
		}
		else
		{
			// (1.0, 0.0)
			vec_ortho[0] = 1.0;
			vec_ortho[1] = 0.0;
		}
	}
	else
	{
		if(fabs(vec[1]) < epsilon)
		{
			// (0.0, 1.0)
			vec_ortho[0] = 0.0;
			vec_ortho[1] = 1.0;
		}
		else 	// most common case
		{
			vec_ortho[0] = 1.0;
			vec_ortho[1] = -vec[0]/vec[1];
		}	
	}

	return vec_ortho;
}

template <>
inline
PointD<3> orthogonalVector<3>(const PointD<3> & vec)
{
	double epsilon = 0.00001;
	PointD<3> vec_ortho(PointD<3>::diagonal(0.0));

	unsigned int zero_coord_count = 0;
	for(std::size_t i = 0; i < 3; i++)
	{
		if(vec[i] < epsilon)
		{
			zero_coord_count++;
		}
	}

	// depending on the number of coordinates equal to zero, there can be trivial cases
	switch(zero_coord_count)
	{
		case 3 :	// zero filled vector, keep vec_ortho as is
			break;
		case 2 :
		{
			vec[0] < epsilon ? vec_ortho[0] = 1.0 : vec_ortho[2] = 1.0;
			break;
		}	
		case 1 :
		{
			// find the zero coordinate
			std::size_t zero_coord_index = 0;
			for(std::size_t i = 0; i < 3; i++)
			{
				if(vec[i] < epsilon)
				{
					zero_coord_index = i;
				}
			}
			vec_ortho[zero_coord_index] = 1.0;
			break;
		}
		case 0 :
		{
			vec_ortho[0] = 0.0;
			vec_ortho[1] = 1.0;
			vec_ortho[2] = -vec[1] / vec[2];
			break;
		}
		default:
		{
			std::cout << "GeomHelpers::orthogonalVector() : DEFAULT SWITCH CASE NOT SUPPOSED TO HAPPEN" << std::endl;
			break;
		}
	}

	return vec_ortho / vec_ortho.norm();
}


}

#endif // !defined GEOMHELPERS_h

#undef GEOMHELPERS_RECURSES
#endif // else defined(GEOMHELPERS_RECURSES)