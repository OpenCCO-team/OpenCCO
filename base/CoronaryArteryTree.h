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
  
  // my_aPerf: Perfusion area
  double my_aPerf =10000;
  
  // End biological parameters
  //-----------------------------
  
  
  
  //-----------------------------
  // Internal algorithm parameter

  //myKTerm: current number of terminal segments
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
    my_aPerf = aPerf;
    my_NTerm = nTerm;
    myDThresold = sqrt(M_PI*myRsupp*myRsupp/myKTerm);
    // Generate the first random terminal point
    Point2D pTerm = ptRoot; //generateRandomPtOnDisk(myTreeCenter, myRsupp);
    
    // Construction of the special root segment
    Segment<Point2D> s;
    s.myRadius = aRadius;
    s.myCoordinate = Point2D(0,myRsupp);//ptRoot;
    s.myLength = 0;
    s.myIndex = 0;
    myVectSegments.push_back(s);
    myVectParent.push_back(0); //if parent index is itsef it is the root (special segment of length 0).
    myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0)); // if children index is itself, it is an end segment.
    my_rPerf = myRsupp;//(ptRoot-myTreeCenter).norm();
    // Construction of the first segment after the root
    Segment<Point2D> s1;
    s1.myRadius = aRadius;
    s1.myCoordinate = pTerm;
    s1.myLength = (ptRoot-pTerm).norm();
    s1.myIndex = 1;
    myVectSegments.push_back(s1);
    myVectTerminals.push_back(1);
    myVectParent.push_back(0); //if parent index is the root
    myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0)); // if children index is itself, it is an end segment.
    DGtal::trace.info() << "Construction initialized..." << std::endl;
  };
  
  
  
  
  // ----------------------- Interface --------------------------------------
  
  bool addFirstSegment(const Point2D &p);
  
  /**
   * Tries to add a new segment from a given point and it parent index.
   * @param p the extremity of the nex segment to be created
   * @param indexParent the index of the parent.
   * @return true of the new segment is created, false in the other case.
   * (for instance if an intersection to previous point was present)
   **/
  
  bool addSegmentFromPoint(const Point2D &p,  unsigned int nearIndex,
                           double rLeft = 1.0, double rRight = 1.0);
  
  
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
   * Computes the n nearest index of neighborhood segments  of a given point
   * @param p : the point considered to get the nearest point
   * @param n : the number of nearest point to be recovered
   */
  
  std::vector<unsigned int> getN_NearestSegments(const Point2D &p,
                                                 unsigned int n) const;
  /**
   * Computes the distance from and the segment represented with the index  and the point given as argument.
   * @param index : the index of the segement used for the comparison

   */
 
  double getDistance(unsigned int index, const Point2D &p ) const;
  
  
  /**
   * Compute if a segment has intersection on the n nearest segments.
   * @param p0 one extremity of one segment
   * @param p1 another extremity of one segment
   * @param n : the number of nearest point to be considered
   */
  
  bool hasNearestIntersections(const Point2D &p0,
                               const Point2D &p1, unsigned int n) const;
  
  
  
  
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
  void kamyiaOptimization(unsigned int index, unsigned int nbIter = 100);
  
  
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






