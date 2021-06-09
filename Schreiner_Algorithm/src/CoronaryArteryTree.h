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

#include "DGtal/helpers/StdDefs.h"

const double pi =   std::atan(1)*4;


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
  template <typename TPoint>
  struct Segment{
    TPoint myCoordinate; // represent the distal point of the segment.
    unsigned int index = 0; // index in the vectSegments
    // radius of the tubular section.
    double radius = 1.0;
    // length of the tubular section.
    double length = 0.0;
    int Generation=0;
  };


  typedef std::pair<unsigned int, unsigned int> SegmentDaugther;
  typedef DGtal::Z2i::RealPoint Point2D;

  std::vector<std::pair<unsigned int, unsigned int> >  myVectDaughters; // to recover the Daugher left (first) and right (second) on an indexed segment.
  std::vector<unsigned int >  myVectParent; // to recover the parent of an indexed segement
  std::vector<Segment<Point2D>> myVectSegments;
  
  
  
  //-----------------------------
  // Global biological parameters
  
  // Number of terminal segments
  unsigned int nTerm = 0;

  // rPerf radius of circular area (m)
  double rPref = 0.05;

  // qPerf LAD flow in vasodilated
  double qPerf = 0.00000833;

  // pPerf perfusion pressure
  double pPerf = 13300.0;

  // pTerm pressure et distal ends of terminal segment
  double pTerm = 4000;
  
  // qTerm flow in one terminal segment
    double qTerm = 1000;
    
  // gamma: bifurcation exponent (2.10-3.0)
    double gamma = 2.10;
    
  // aPerf: Perfusion area
    double aPerf =10000;
  
  // End biological parameters
  //-----------------------------


  
  //-----------------------------
  // Internal algorithm parameter

  unsigned int myKTerm = 1;
  
  double myDThresold = 0;
  
  double myRsupp = 0;
    
  unsigned int Nterm=1;
    
    unsigned int MaxGene=1;
  
  Point2D myTreeCenter;
  
  
  
  // End: Internal algorithm parameter
  //-----------------------------
  
  
  
public: 
 
  /**
   * Default constructor
   **/
  
  CoronaryArteryTree(const Point2D &ptRef, double aPref, unsigned int Nterm){
    myRsupp = sqrt((aPref/Nterm)/pi);
    aPerf = aPref;
    Nterm=Nterm;
    myDThresold=sqrt(pi*myRsupp*myRsupp/myKTerm);
    
    myTreeCenter = Point2D(0,0);
    Segment<Point2D> s;
    s.radius=15.0; s.myCoordinate = ptRef; s.length=5.0;
    s.index = 0; 
    myVectSegments.push_back(s);
    myVectParent.push_back(0); //if parent index is itsef it is the root (special segment of length 0).
    myVectDaughters.push_back(std::pair<unsigned int, unsigned int>(0,0)); // if daugthers index is itself, it is an end segment.
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
  
  bool addSegmentFromPoint(const Point2D &p,  unsigned int nearIndex);
  
   
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
    
  void exportDisplay(const std::string &fileName = "result.eps");

  Point2D getSegmentCenter(const Segment<Point2D> &s);

  Point2D getSegmentCenter(unsigned int i);
    double GetTotalVolume(const Point2D &p1,const Point2D &p2,const Point2D &p3,const Point2D &pOpti);
  unsigned int getNearestSegment(const Point2D &pt);
  
  unsigned int getParentSegment(const Segment<Point2D> &s);
  unsigned int getDaughterLeft(const Segment<Point2D> &s);
  unsigned int getDaughterRigth(const Segment<Point2D> &s);
    bool addSegment(const Point2D &NewPoint,const Point2D &OptimizePoint, unsigned int nearIndex);
  /**
   * Computes the distance criteria (d_crit) computed from the orthogonal projection or distance to end point.
   *
   */
  double compDistCriteria(const Point2D &p, unsigned int indexNode);
  double dProjCalculation(const Point2D &p,unsigned int Index );
  double dCritCalculation(const Point2D &p,unsigned int Index );
  Point2D FindBarycenter(const Point2D &p, unsigned int index);
  bool updateRsupp();
  bool updateTreshold();
  double GetLength(unsigned int Index);
  bool updateRadius();
    bool updateRadius2(unsigned int index );
    Point2D FindOptimalPosition(unsigned int Index,const Point2D &p);
    bool updateScale( double resize_factor);
    bool updateGeneration();
    int FindOptimalSegment(const Point2D &p);
    double FindXmax(int xDim, int yDim);
    double FindYmax(int xDim, int yDim);
    int AddFirstSegmentonImage();
//    double FindEdge(int xDim, int yDim);
    Point2D fromCircleToImage(std::string fileName, double x, double y, int xdim,int ydim );
    Point2D fromImageToCircle(int ximage,int yimage,int xdim,int ydim);

    bool hasIntersections(Segment<CoronaryArteryTree::Point2D> S1, Point2D newPoint);
  /**
   * Use to display object.
   * @param out  the object of class 'GeodesicGraphComputer' to write.
   */
  void selfDisplay( std::ostream & out ) const;

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


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CORONARY_ARTERY_TREE_h

#undef CORONARY_ARTERY_TREE_RECURSES
#endif // else defined(CORONARY_ARTERY_TREE_RECURSES)






