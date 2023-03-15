
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
        for(auto i=0; i<TPoint::dimension; i++){ ptProjected[i]=ptA[i];}
        return true;
    }
    if (ptB==ptC)
    {
        for(auto i=0; i<TPoint::dimension; i++){ ptProjected[i]=ptB[i];}
        return true ;
    }
    
    TPointD vAB = ptB - ptA;//(ptB[0]- ptA[0], ptB[1]- ptA[1]);
    TPointD vABn = vAB / vAB.norm();//((double)vAB[0], (double)vAB[1]);
    //double norm = vABn.norm();
    //vABn[0] /= norm;
    //vABn[1] /= norm;
    
    TPointD vAC = ptC-ptA;// (ptC[0]-ptA[0], ptC[1]-ptA[1]);
    double distPtA_Proj = vAC.dot(vABn);
    
    for(auto i=0; i<TPoint::dimension; i++){ ptProjected[i]= ptA[i]+vABn[i]*(distPtA_Proj);}
    //ptProjected[0] = ptA[0]+vABn[0]*(distPtA_Proj);
    //ptProjected[1] = ptA[1]+vABn[1]*(distPtA_Proj);
    bool res = false;
    for(auto i=0; i<TPoint::dimension; i++) { res = res || (ptA[i]<ptB[i] && ptProjected[i]<=ptB[i]); }
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
segment2segmentDistance(const DGtal::PointVector<TDim, double> &segA, const DGtal::PointVector<TDim, double> &segB,
                const DGtal::PointVector<TDim, double> &segC, const DGtal::PointVector<TDim, double> &segD)
{
  DGtal::PointVector<TDim, double> A, B, C, D;
  DGtal::PointVector<TDim, double> d1, d2, d12;
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
point2segmentDistance(const DGtal::PointVector<TDim, double> &segA, const DGtal::PointVector<TDim, double> &segB,
                const DGtal::PointVector<TDim, double> &ptP)
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
isIntersecting(const DGtal::PointVector<TDim, double> &segA, const DGtal::PointVector<TDim, double> &segB, double rAB,
                      const DGtal::PointVector<TDim, double> &segC, const DGtal::PointVector<TDim, double> &segD, double rCD)
{
  DGtal::PointVector<TDim, double> cAB = (segA + segB)/2.0;
  DGtal::PointVector<TDim, double> cCD = (segC + segD)/2.0;
  double lAB = (segA - segB).norm();
  double lCD = (segC - segD).norm();
  double dC = (cAB - cCD).norm();
  bool res = dC > ((lAB + lCD)/2.0 + rAB + rCD);
  if(res==true)
    return false;
  double distanceSeg = segment2segmentDistance<TDim>(segA, segB, segC, segD);
  if(distanceSeg > (rAB + rCD))
    return false;
  return true;
}




template<typename TPoint, typename TDSet>
inline
TDSet
pointsOnSphere(const TPoint & ptCenter,
               double radius)
{
    typedef  DGtal::SpaceND<TPoint::dimension,
                            typename TPoint::Component> Space;
    typedef DGtal::HyperRectDomain<Space> Dom;
    
    typedef DGtal::ImplicitBall< Space > MyBall;
    MyBall disk( ptCenter, radius-0.5 );
    MyBall diskDilate( ptCenter, radius+0.5 );
    typedef DGtal::EuclideanShapesCSG< MyBall, MyBall > Minus;
    Minus border ( diskDilate );
    border.minus( disk );
    typedef DGtal::GaussDigitizer< Space, Minus > MyGaussDigitizer;
    MyGaussDigitizer digShape;
    digShape.attach( border );
    digShape.init( border.getLowerBound(), border.getUpperBound(), 1 );
    Dom domainShape = digShape.getDomain();
    TDSet  aSet(domainShape);
    DGtal::Shapes<Dom>::digitalShaper( aSet, digShape );
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
    typedef typename TImage::Domain::Space Space;
    typedef DGtal::ExactPredicateLpSeparableMetric<Space, TImage::Domain::dimension> L2Metric;
    
    TImageDistance res (image.domain());
    typedef DGtal::functors::IntervalForegroundPredicate<TImage> Binarizer;
    typedef DGtal::DistanceTransformation<Space, Binarizer, L2Metric> DTL;
    L2Metric l2metric;
    Binarizer b(image, threshold, 255);
    DTL dt(&image.domain(),&b, &l2metric);
    for (auto p: dt.domain()){
        res.setValue(p, dt(p));
    }
    return res;
}


}



#endif // !defined GEOMHELPERS_h

#undef GEOMHELPERS_RECURSES
#endif // else defined(GEOMHELPERS_RECURSES)
