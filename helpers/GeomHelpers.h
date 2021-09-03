


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


namespace GeomHelpers {

//template<typename TPoint>
//inline
//TPoint
//generateRandomPtOnDisk(const TPoint &ptCenter, double r)
//{
//  double a = (rand()%360)/(2.0*M_PI);
//  double rR = ((double)rand() / RAND_MAX)*r;
//  return TPoint(ptCenter[0]+rR*cos(a), ptCenter[1]+rR*sin(a));
//}
template<typename TPoint>
inline
TPoint
generateRandomPtOnDisk(const TPoint &ptCenter, double r)
{
  bool found = false;
  double x = 0.0;
  double y = 0.0;
  
  while(!found){
    x =  ((double)rand() / RAND_MAX)*2.0*r - r;
    y =  ((double)rand() / RAND_MAX)*2.0*r - r;
    found = x*x + y*y < r*r;
  }
  return TPoint(x+ptCenter[0], y+ptCenter[1]);
}



template<typename TPoint, typename TImage, typename TImageDist>
inline
TPoint
generateRandomPtOnImageDomain(const TImage &image, unsigned int fgTh,
                              const TImageDist &imageDistance,
                              unsigned int nbTry = 100)
{
  bool found = false;
  unsigned int x = 0;
  unsigned int y = 0;
  TPoint pMin = image.domain().lowerBound();
  TPoint pMax = image.domain().upperBound();
  int dx = pMax[0] - pMin[0];
  int dy = pMax[1] - pMin[1];
  DGtal::Z2i::Point pCand;
  unsigned int n = 0;
  while(!found && n < nbTry){
    x =  rand()%dx;
    y =  rand()%dy;
    pCand[0] = pMin[0] +x;
    pCand[1] = pMin[1] +y;
    found = image(pCand)>=fgTh && abs(imageDistance(pCand)) >= 10.0;
    n++;
  }
  if (n >= nbTry){
    for(auto p : image.domain()){if (image(p)>=fgTh && abs(imageDistance(p)) >= 10.0 ) return p;}
  }
  return pCand;
}



/**
 * Check if the segment defined by two points intersect the domain.
 *
 * @param image defining the domain
 * @param fgTh the threshol defining the foreground
 * @param pt1 first point of the segment
 * @param pt2  second point of the segment
 */
template<typename TImage>
inline
bool
checkNoIntersectDomain(const TImage &image, unsigned int fgTh,
                       const DGtal::Z2i::Point &pt1,
                       const DGtal::Z2i::Point &pt2)
{
  if ( !image.domain().isInside(pt1) ||
      !image.domain().isInside(pt2)){
    return false;
  }
  DGtal::Z2i::RealPoint dir = pt2 - pt1;
  dir /= dir.norm();
  DGtal::Z2i::RealPoint p (pt1[0], pt1[1]);
  for (unsigned int i = 0; i<(pt2 - pt1).norm(); i++){
    DGtal::Z2i::RealPoint p = pt1+dir*i;
    if (image(DGtal::Z2i::Point(static_cast<int>(p[0]),static_cast<int>(p[1]) )) < fgTh)
      return false;
  }
  
  return true;
}

/**
 * Dertermines if a point is on the right of a line represented by two points [ptA, ptB]
 */
template<typename TPoint>
inline
bool
isOnRight(const TPoint &ptA, const TPoint &ptB, const TPoint &ptC ){
  auto u = ptB-ptA;
  auto v = ptC-ptA;
  return u[0]*v[1]- u[1]*v[0] < 0.0;
}


template<typename TPoint>
inline
bool
isInsideCircle(const TPoint &ptCenter,const TPoint &p,  double radius){
  return ((p[0]-ptCenter[0])*(p[0]-ptCenter[0])+(p[1]-ptCenter[1])*(p[1]-ptCenter[1])) <= radius*radius;
}

template<typename TPoint>
inline
bool
isInsideSphere(const TPoint &ptCenter,const TPoint &p,  double radius){
  return ((p[0]-ptCenter[0])*(p[0]-ptCenter[0])+(p[1]-ptCenter[1])*(p[1]-ptCenter[1])+(p[2]-ptCenter[2])*(p[2]-ptCenter[2])) <= radius*radius;
}



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

template<typename TPoint, typename TPointD>
inline
bool
projectOnStraightLine(const TPoint & ptA,
                      const TPoint & ptB,
                      const TPoint & ptC,
                      TPointD & ptProjected)
{
  if (ptA==ptC)
  {
    ptProjected=ptA;
    return true;
  }
  if (ptB==ptC)
  {
    ptProjected=ptB;
    return true ;
  }
  
  TPointD vAB (ptB[0]- ptA[0], ptB[1]- ptA[1]);
  TPointD vABn ((double)vAB[0], (double)vAB[1]);
  double norm = vABn.norm();
  vABn[0] /= norm;
  vABn[1] /= norm;
  
  TPointD vAC (ptC[0]-ptA[0], ptC[1]-ptA[1]);
  double distPtA_Proj = vAC.dot(vABn);
  
  ptProjected[0]= ptA[0]+vABn[0]*(distPtA_Proj);
  ptProjected[1] = ptA[1]+vABn[1]*(distPtA_Proj);
  
  return  distPtA_Proj>=0 && ((ptA[0]<ptB[0] && ptProjected[0]<=ptB[0] ) ||
                              (ptA[0]>ptB[0] && ptProjected[0]>=ptB[0] ) ||
                              (ptA[0]==ptB[0] && ptA[1]<ptB[1] && ptProjected[1]<=ptB[1]) ||
                              (ptA[0]==ptB[0] && ptA[1]>=ptB[1] && ptProjected[1]>=ptB[1]));
}


template<typename TPoint>
inline
DGtal::Z2i::DigitalSet
pointsOnCircle(const TPoint & ptCenter,
               double radius)
{
  typedef DGtal::ImplicitBall< DGtal::Z2i::Space > MyBall;
  MyBall disk( ptCenter, radius-0.5 );
  MyBall diskDilate( ptCenter, radius+0.5 );
  typedef DGtal::EuclideanShapesCSG< MyBall, MyBall > Minus;
  Minus border ( diskDilate );
  border.minus( disk );
  typedef DGtal::GaussDigitizer< DGtal::Z2i::Space, Minus > MyGaussDigitizer;
  MyGaussDigitizer digShape;
  digShape.attach( border );
  digShape.init( border.getLowerBound(), border.getUpperBound(), 1 );
  DGtal::Z2i::Domain domainShape = digShape.getDomain();
  DGtal::Z2i::DigitalSet aSet( domainShape );
  DGtal::Shapes<DGtal::Z2i::Domain>::digitalShaper( aSet, digShape );
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
  if ( d==0.0 )
  {
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
static bool kamyiaOpt(double gamma, double deltaP1, double deltaP2, double f0, double f1, double f2, double l0, double l1, double l2, double &xx1, double &xx2) {
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
  CostFunction* cost_function =
  new AutoDiffCostFunction<CostOptPos, 2, 1, 1>(f);
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
getImageDistance(const TImage &image, unsigned int threshold=128){
  TImageDistance res (image.domain());
  typedef DGtal::functors::IntervalForegroundPredicate<TImage> Binarizer;
  typedef DGtal::DistanceTransformation<DGtal::Z2i::Space, Binarizer, DGtal::Z2i::L2Metric> DTL2;
  Binarizer b(image, threshold, 255);
  DTL2 dt(&image.domain(),&b, &DGtal::Z2i::l2Metric);
  for (auto p: dt.domain()){
    res.setValue(p, dt(p));
  }
  return res;
}
}


#endif // !defined GEOMHELPERS_h

#undef GEOMHELPERS_RECURSES
#endif // else defined(GEOMHELPERS_RECURSES)
