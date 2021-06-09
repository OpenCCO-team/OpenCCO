#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"



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
  vABn = vABn/vABn.norm();
  TPointD vAC (ptC[0]-ptA[0], ptC[1]-ptA[1]);
  double distPtA_Proj = vAC.dot(vABn);

  ptProjected[0]= ptA[0]+vABn[0]*(distPtA_Proj);
  ptProjected[1] = ptA[1]+vABn[1]*(distPtA_Proj);

  return  distPtA_Proj>=0 && ((ptA[0]<ptB[0] && ptProjected[0]<=ptB[0] ) ||
                              (ptA[0]>ptB[0] && ptProjected[0]>=ptB[0] ) ||
                              (ptA[0]==ptB[0] && ptA[1]<ptB[1] && ptProjected[1]<=ptB[1]) ||
                              (ptA[0]==ptB[0] && ptA[1]>=ptB[1] && ptProjected[1]>=ptB[1]));
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





