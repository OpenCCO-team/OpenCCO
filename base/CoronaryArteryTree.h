#if defined(CORONARY_ARTERY_TREE_RECURSES)
#error Recursive header files inclusion detected in CoronaryArteryTree.h
#else // defined(CORONARY_ARTERY_TREE_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CORONARY_ARTERY_TREE_RECURSES

#if !defined CORONARY_ARTERY_TREE_h
/** Prevents repeated inclusion of headers. */
#define CORONARY_ARTERY_TREE_h



#include <iostream>
#include <math.h>
#include "DGtal/io/boards/Board2D.h"

#include "DGtal/helpers/StdDefs.h"
#include "geomhelpers.h"


#include "ceres/ceres.h"
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;


///////////////////////////////////////////////////////////////////////////////
// class CoronaryArteryTree
/**
 * Description of class 'CoronaryArteryTree' <p>
 *
 * @brief Class to represent and construct a Coronary Artery Tree using throught the CCO algorithm.
 *
 *
 */



class CoronaryArteryTree{
  
  /**
   * Class to handle the reprensetation of the coronary tree.
   */
public:
  
  // Represent the left and right
  typedef std::pair<unsigned int, unsigned int> SegmentChildren;
  typedef DGtal::Z2i::RealPoint Point2D;
  
  template <typename TPoint>
  struct Segment{
    TPoint myCoordinate; // represents the distal point of the segment.
    unsigned int myIndex = 0; // index in the vectSegments
    // radius of the tubular section.
    double myRadius = 1.0;
    // length of the tubular section.
    double myLength = 0.0;
    
    // number of terminal segments in its children segment
    unsigned int myKTerm = 1;
    // hydrodynamc registance (R star)
    double myHydroResistance = 0.00;
    // flow (Qi)
    double myFlow = 0.00;
    // radius ratio of the segment and his brother (ri/rj)
    double myRaidusRatio = 1.0;
    // radius ratio of the segment and its parent
    double beta = 1.0;
    
  };
  
  // to recover the Children left (first) and right (second) on an indexed segment.
  std::vector< SegmentChildren >  myVectChildren;
  
  // to recover the parent of an indexed segement
  std::vector<unsigned int >  myVectParent;
  // represents all the vertices of the graph
  std::vector<Segment<Point2D>> myVectSegments;
  
  
  //-----------------------------
  // Global biological parameters
  
  // my_NTerm: number of terminal segments
  unsigned int my_NTerm = 1;
  
  // my_rPref radius of circular area (m)
  double my_rPerf = 1.0;
  
  // my_qPerf LAD flow in vasodilated
  double my_qPerf = 0.00000833;
  
  // my_pPerf perfusion pressure
  double my_pPerf = 13300.0;
  
  // my_pTerm pressure et distal ends of terminal segment
  double my_pTerm = 4000;
  
  // my_qTerm flow in one terminal segment
  double my_qTerm = 1000;
  
  // my_gamma: bifurcation exponent (2.10-3.0)
  double my_gamma = 3.0;
  
  // my_mu: viscosity of blood
  double my_mu = 3.6;
  
  // my_aPerf: Perfusion area
  double my_aPerf =10000;
  
  // End biological parameters
  //-----------------------------
  
  
  
  //-----------------------------
  // Internal algorithm parameter

  //myKTerm: current number of terminal segments of the tree
  unsigned int myKTerm = 1;
  
  //myDThresold: threshold on the distance criterion of adding a segment
  double myDThresold = 0.0;
  
  //myRsupp: average radius of blackboxes
  double myRsupp = 0.0;
  
  // myTreeCenter: coordinate of the tree center used to define the main domain.
  Point2D myTreeCenter;
  
  
  // myCurrAPerf : represents the current perfusion surface
  double myCurrAPerf = 1.0;
  
  // myVectTerminals : store the index of the terminal segments
  std::vector<unsigned int> myVectTerminals;
  
  // End: Internal algorithm parameter
  //-----------------------------
  
  // To handle display
  DGtal::Board2D myBoard;

  
public: 
  
  /**
   * Default constructor.
   * It generates the first root segment by randomly choose the first terminal point.
   * @param ptRoot: coordinates of the root special segment (who have no parent)
   * @param aPerf: surface of the perfusion.
   * @param nTerm: number of terminal segments.
   **/
  
  CoronaryArteryTree(const Point2D &ptRoot, double aPerf, unsigned int nTerm,
                     double aRadius = 1.0 ){
     
    myTreeCenter = Point2D(0,0);
    myRsupp = sqrt((aPerf/nTerm)/M_PI);
    my_rPerf = myRsupp;
    my_aPerf = aPerf;
    my_NTerm = nTerm;
    myDThresold = sqrt(M_PI*myRsupp*myRsupp/myKTerm);
    // Generate the first random terminal point
    //Point2D pTerm = ptRoot; //generateRandomPtOnDisk(myTreeCenter, myRsupp);
    Point2D pTerm = ptRoot;
    
    // Construction of the special root segment
    Point2D ptRootNew = Point2D(0,myRsupp);
    Segment<Point2D> s;
    s.myRadius = aRadius;
    s.myCoordinate = ptRootNew;
    s.myLength = 0;
    s.myIndex = 0;
    s.myKTerm = 0;
    myVectSegments.push_back(s);
    myVectParent.push_back(0); //if parent index is itsef it is the root (special segment of length 0).
    myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0)); // if children index is itself, it is an end segment.
    
    // Construction of the first segment after the root
    Segment<Point2D> s1;
    s1.myRadius = aRadius;
    s1.myCoordinate = pTerm; //myTreeCenter
    s1.myLength = (ptRootNew-pTerm).norm(); //(ptRootNew-myTreeCenter).norm();
    s1.myIndex = 1;
    s1.myKTerm = 1; //it contains terminal itself
    s1.myHydroResistance = 8.0*my_mu*s1.myLength/M_PI;
    s1.myFlow = my_qTerm;
    s1.myRaidusRatio = 0.0;
    s1.beta = 1.0;
    
    myVectSegments.push_back(s1);
    myVectTerminals.push_back(1);
    myVectParent.push_back(0); //if parent index is the root
    myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0)); // if children index is itself, it is an end segment.
    DGtal::trace.info() << "Construction initialized..." << std::endl;
  };
  
  
  
  
  // ----------------------- Interface --------------------------------------
  
  bool addFirstSegment(const Point2D &p);
  
  double computeTotalVolume(unsigned int segIndex);
  
  bool isAddable(const Point2D &p, unsigned int segIndex, unsigned int nbIter);
  
  bool isIntersecting(const Point2D &p,  unsigned int nearIndex, unsigned int nbNeibour = 10, double minDistance = 5.0);
  /**
   * Tries to add a new segment from a given point and it parent index.
   * @param p the extremity of the nex segment to be created
   * @param nearIndex the index of the near segement to p.
   * @return true of the new segment is created, false in the other case.
   * (for instance if an intersection to previous point was present)
   **/
  
  bool addSegmentFromPoint(const Point2D &p,  unsigned int nearIndex,
                           double rLeft = 1.0, double rRight = 1.0);
  
  bool addSegmentFromPointBK(const Point2D &p,  unsigned int nearIndex,
                           double rLeft = 1.0, double rRight = 1.0);
  
  
  /**
   * Update the distribution of segmental flows after adding a new segment (new bifurcation)
   * @param segIndex index of the parent segment to be updated
   */
  void updateFlowTerminal(unsigned int segIndex);
  void updateFlowParameters(unsigned int segIndex);
  
  /**
   * Update ...
   * @param segIndex index of the parent segment to be updated
   */
  void updateRadius(unsigned int segIndex);
  
  /**
   * Update the root radius after updating flow parameters (as given in Eq. 18)
   */
  void updateRootRadius();
  
  /**
   * Update the segment radius after updating flow parameters (as given in Eq. 19)
   */
  void updateSegmentRadiusToRoot(unsigned int segIndex);
  
  /**
   * Compute the total tree volume (as given in Eq. 20)
   */
  double computeTreeVolume(double mu, double lambda);
  
  /**
   * Tries to add a new segment from a given point.
   * @param p the extremity of the nex segment to be created
   * @param indexParent the index of the parent.
   * @return true of the new segment is created, false in the other case.
   * (for instance if an intersection to previous point was present)
   **/
  
  bool addSegmentFromPoint(const Point2D &p);
  bool addSegmentFromPointWithBarycenter(const Point2D &p);
  bool addSegmentFromPointWithBarycenter(const Point2D &p, unsigned int nearIndex);
  /**
   * Export the current display of the tree.
   */
  
  void boardDisplay(bool clearDisplay = true);
  void exportBoardDisplay(const std::string &fileName = "result.eps",
                          bool updateDisplay = true );
  
  
  Point2D getSegmentCenter(const Segment<Point2D> &s);
  
  Point2D getSegmentCenter(unsigned int i);
  double GetTotalVolume(const Point2D &p1,const Point2D &p2,const Point2D &p3,const Point2D &pOpti);
  unsigned int getNearestSegment(const Point2D &pt);
  
  unsigned int getParentSegment(const Segment<Point2D> &s);
  unsigned int getLeftChild(const Segment<Point2D> &s);
  unsigned int getRightChild(const Segment<Point2D> &s);
  bool addSegment(const Point2D &NewPoint,const Point2D &OptimizePoint, unsigned int nearIndex);
  /**
   * Computes the distance criteria (d_crit) computed from the orthogonal projection or distance to end point.
   *
   */
  double compDistCriteria(const Point2D &p, unsigned int indexNode);
  double dProjCalculation(const Point2D &p,unsigned int Index );
  double dCritCalculation(const Point2D &p,unsigned int Index );
  Point2D FindBarycenter(const Point2D &p, unsigned int index);
  void updateTreshold();
  double GetLength(unsigned int Index);
  bool updateRadius();
  bool updateRadius2(unsigned int index );
  double FindXmax(int xDim, int yDim);
  double FindYmax(int xDim, int yDim);
  int AddFirstSegmentonImage();
  Point2D fromCircleToImage(std::string fileName, double x, double y, int xdim,int ydim );
  Point2D fromImageToCircle(int ximage,int yimage,int xdim,int ydim);
  
  /// New from updating code... (BK+PN)
  /**
   * From a segment returns a vector of segment index representing the path to the root.
   */
  std::vector<unsigned int> getPathToRoot(const Segment<Point2D> &s);
  
  void udpatePerfusionArea();
  
  /**
   * Generate a new location with distance constraints.
   * @param nbTrials: number of trials before reducing the distance constaint value
   *
   */
  Point2D generateNewLocation(unsigned int nbTrials = 1000);

  std::pair<Point2D, bool> generateALocation();

  /**
   * Computes the distance from a segment represented with the index  and the point given as argument.
   * @param index : the index of the segement used for the comparison
   * @param p : a point
   */
  double getDistance(unsigned int index, const Point2D &p ) const;
  
  
  
  /**
   * Computes the projected distance from a segment represented with the index  and the point given as argument.
   * @param index : the index of the segement used for the comparison
   * @param p : a point
   */
  double getProjDistance(unsigned int index, const Point2D &p) const;
  /**
   * Computes the projected distance from a segment represented with the index  and the point given as argument.
   * @param p0 : a point representing one extremity
   * @param p1 : a point representing another extremity
   * @param p : a point to be projected
   */
  double getProjDistance(const Point2D &p0, const Point2D &p1, const Point2D &p) const;
  
  /**
   * Check if a new added point is too close the nearest segment.
   * @param p : a point
   * @param minDist: min distance
   * @return true is a point is too close
   */
  bool isToCloseFromNearest(const Point2D &p, double minDist) const;
  
  
  
  /**
   * Computes the n nearest index of neighborhood segments  of a given point
   * @param p : the point considered to get the nearest point
   * @param n : the number of nearest point to be recovered
   */
  std::vector<unsigned int> getN_NearestSegments(const Point2D &p,
                                                 unsigned int n) const;
  
  /**
   * Compute if a segment has intersection on the n nearest segments.
   * @param p0 one extremity of one segment
   * @param p1 another extremity of one segment
   * @param n : the number of nearest point to be considered
   */
  
  bool hasNearestIntersections(const Point2D &p0,
                               const Point2D &p1, unsigned int n) const;
  
  
  
  /**
   * Computes if an bifurcation has intersections on the n nearest segments.
   * It uses the middle point
   * @param indexPFather index extremity of initial segment (father).
   * @param indexPChild index  extremity  of initial segment (child).
   * @param pAdded point added at the origin of the creation of bifurcation.
   * @param pBifurcation the central point of the bifurcation.
   * @param n : the number of nearest point to be considered
   */
  
  bool hasNearestIntersections(unsigned int indexPFather,
                               unsigned int indexPChild,
                               const Point2D &pAdded,
                               const Point2D &pBifurcation, unsigned int n) const;
  
  
  
  
  /// Fin New from updating code... (BK+PN)
  
  
  
  
  bool hasIntersections(Segment<Point2D> S1, Point2D newPoint);
  /**
   * Use to display object.
   * @param out  the object of class 'GeodesicGraphComputer' to write.
   */
  void selfDisplay( std::ostream & out ) const;
  
  
  /**
   * @param index: the segment index
   */
  bool kamyiaOptimization(unsigned int index, unsigned int nbIter = 100);
  
  bool kamyiaOptimization(const DGtal::Z2i::RealPoint& pParent, const Segment<Point2D>& sCurrent, const Segment<Point2D>& sL, const Segment<Point2D>& sR, unsigned int nbIter, DGtal::Z2i::RealPoint& pOpt, double& r0, double& r1, double& r2);
  
};



/**
 * Overloads 'operator<<' for displaying objects of class 'CoronaryArteryTree'.
 * @param out the output stream where the object is written.
 * @param aCoronaryTree the object of class 'CoronaryArteryTree' to write.
 * @return the output stream after the writing.
 */

std::ostream&
operator<< ( std::ostream & out, const CoronaryArteryTree & aCoronaryTree );

bool
operator==(CoronaryArteryTree::Segment<CoronaryArteryTree::Point2D>  S1, CoronaryArteryTree::Segment<CoronaryArteryTree::Point2D> S2);


bool
operator!=(CoronaryArteryTree::Segment<CoronaryArteryTree::Point2D>  S1, CoronaryArteryTree::Segment<CoronaryArteryTree::Point2D> S2);







///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods 


///////////////////////////////////////////////////////////////////////////////

#endif // !defined CORONARY_ARTERY_TREE_h

#undef CORONARY_ARTERY_TREE_RECURSES
#endif // else defined(CORONARY_ARTERY_TREE_RECURSES)






