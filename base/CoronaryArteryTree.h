#pragma once

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
#include "DGtal/images/ImageSelector.h"

#include "DGtal/helpers/StdDefs.h"

#include "GeomHelpers.h"


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

template <int TDim>
class CoronaryArteryTree{
  
  /**
   * Class to handle the reprensetation of the coronary tree.
   */
public:
  
  
  // Domain
  typedef DGtal::SpaceND< TDim, int >   SpaceCT;
  typedef DGtal::HyperRectDomain<SpaceCT> DomCT;
  typedef DGtal::PointVector<TDim, int> TPoint;
  typedef DGtal::PointVector<TDim, double> TPointD;
  
  // Represent the left and right
  typedef std::pair<unsigned int, unsigned int> SegmentChildren;

  typedef typename DGtal::ImageSelector < DomCT, unsigned char>::Type Image;
  typedef DGtal::ImageContainerBySTLVector< DomCT, int> ImageDist;

  struct Segment{
    // Distal point of the segment.
    TPointD myCoordinate;
    // index in the vectSegments.
    unsigned int myIndex = 0;
    
    // radius of the tubular section.
    double myRadius = 1.0;
    // number of terminal segments in its children segment
    unsigned int myKTerm = 1;
    // hydrodynamc registance (R star)
    double myResistance = 0.00;
    // flow (Qi)
    double myFlow = 0.00;
    // radius ratio of the segment and its parent
    double myBeta = 1.0;
  };
  
  // to recover the Children left (first) and right (second) on an indexed segment.
  std::vector< SegmentChildren >  myVectChildren;
  
  // to recover the parent of an indexed segement
  std::vector<unsigned int >  myVectParent;
  // represents all the vertices of the graph
  std::vector< Segment > myVectSegments;
  // to store the index of the terminal segments
  std::vector<unsigned int> myVectTerminals;
  
  //-----------------------------
  // Global biological parameters
  
  // my_NTerm: number of terminal segments
  unsigned int my_NTerm = 1;
  
  // my_rPref final radius of circular area (m)
  double my_rPerf = 1.0;
  
  // my_qPerf LAD flow in vasodilated
  double my_qPerf = 8.33e3; //0.00000833;
  
  // my_pPerf perfusion pressure
  double my_pPerf = 1.33e4; //13300.0;
  
  // my_pTerm pressure et distal ends of terminal segment
  double my_pTerm = 8.40e3;
  
  // my_qTerm flow in one terminal segment
  double my_qTerm = 1000;
  
  // my_gamma: bifurcation exponent (2.10-3.0)
  double my_gamma = 3.0;
  
  // my_mu: viscosity of blood
  double my_mu = 3.6e-3; //3.6;
  
  // my_aPerf: Perfusion area
  double my_aPerf = 10000;
  
  // my_pDrop = my_pPerf-my_pTerm
  double my_pDrop = 4900;
  
  // End biological parameters
  //-----------------------------
  
  //-----------------------------
  // Internal algorithm parameter

  //myKTerm: current number of terminal segments of the tree
  unsigned int myKTerm = 1;
  
  //myRsupp: average radius of blackboxes
  double myRsupp = 0.0;
  
  // myTreeCenter: coordinate of the tree center used to define the main domain.
  TPointD myTreeCenter;
  
  // myNumNeighbor : represents the number of nearest neighbours to be tested
  int myNumNeighbor = 20;
  
  // myLengthFactor : scale factor (updated during tree growth after each added bifurcation)
  double myLengthFactor = 1.0;
  
  // End: Internal algorithm parameter
  //-----------------------------

  
  
  
  // To handle display
  DGtal::Board2D myBoard;

  
protected:
  Image myImageDomain = Image(DomCT());
  ImageDist myImageDist = ImageDist(DomCT());
  unsigned int myForegroundThreshold = 128;
  bool myIsImageDomainRestrained = false;

  
  
public: 
  
  /**
   * @brief Default constructor.
   * @brief It generates the first root segment with tree center as the first terminal point.
   * @param aPerf: surface of the perfusion.
   * @param nTerm: number of terminal segments.
   **/
  
  CoronaryArteryTree(double aPerf, unsigned int nTerm, double aRadius = 1.0,
                     TPointD treeCenter = TPointD::diagonal(0) ){
    assert(nTerm>=1);
    for (auto i=0; i < TDim; i++){myTreeCenter[i]=treeCenter[i];}
    if(TDim==2) {
      myRsupp = sqrt(aPerf/(nTerm*M_PI));
      my_rPerf = sqrt(aPerf/M_PI);
    }
    else {//TDim==3
      myRsupp = pow(3.0*aPerf/(4.0*M_PI*nTerm),1.0/3.0);
      my_rPerf = pow(3.0*aPerf/(4.0*M_PI),1.0/3.0);
    }
    my_aPerf = aPerf;
    my_NTerm = nTerm;
    my_qTerm = my_qPerf / my_NTerm;
    if(nTerm > 250)
      my_pTerm = 7.98e3;
    my_pDrop = my_pPerf-my_pTerm;
    //myDThresold = sqrt(M_PI*myRsupp*myRsupp/myKTerm);
    
    // Construction of the special root segment
    TPointD ptRoot = myTreeCenter;
    ptRoot[1] += my_rPerf;// Point2D(0,my_rPerf);
    Segment s;
    s.myRadius = aRadius;
    s.myCoordinate = ptRoot;
    s.myIndex = 0;
    s.myKTerm = 0;
    myVectSegments.push_back(s);
    myVectParent.push_back(0); //if parent index is itsef it is the root (special segment of length 0).
    myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0)); // if children index is itself, it is an end segment.
    updateLengthFactor();
    
    // Construction of the first segment after the root
    Segment s1;
    s1.myRadius = aRadius;
    s1.myCoordinate = myTreeCenter;
    double myLength = (ptRoot-s1.myCoordinate).norm()*myLengthFactor;
    s1.myIndex = 1;
    s1.myKTerm = 1; //it contains terminal itself
    s1.myResistance = 8.0*my_mu*myLength/M_PI;
    s1.myFlow = my_qTerm;
    s1.myBeta = 1.0;
    
    myVectSegments.push_back(s1);
    myVectTerminals.push_back(1);
    myVectParent.push_back(0); //if parent index is the root
    myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0)); // if children index is itself, it is an end segment.
    updateRootRadius();
    DGtal::trace.info() << "Construction initialized..." << std::endl;
  };
  
  /**
   * @brief Constructor.
   * @brief It generates the first root segment by randomly choose the first terminal point.
   * @param ptRoot: coordinates of the root special segment (who have no parent)
   * @param aPerf: surface of the perfusion.
   * @param nTerm: number of terminal segments.
   **/
  /*
  CoronaryArteryTree(const TPointD &ptRoot, double aPerf, unsigned int nTerm, double aRadius = 1.0 ){
    assert(nTerm>=1);
    myTreeCenter = TPointD::diagonal(0.0);
    if(TDim==2) {
      myRsupp = sqrt(aPerf/(nTerm*M_PI));
      my_rPerf = sqrt(aPerf/M_PI);
    }
    else {//TDim==3
      myRsupp = pow(3.0*aPerf/(4.0*M_PI*nTerm),1.0/3.0);
      my_rPerf = pow(3.0*aPerf/(4.0*M_PI),1.0/3.0);
    }
    my_aPerf = aPerf;
    my_NTerm = nTerm;
    my_qTerm = my_qPerf / my_NTerm;
    if(nTerm > 250)
      my_pTerm = 7.98e3;
    my_pDrop = my_pPerf-my_pTerm;
    //myDThresold = sqrt(M_PI*myRsupp*myRsupp/myKTerm);
    
    // Construction of the special root segment
    assert((ptRoot - myTreeCenter).norm() == my_rPerf); //ptRoot must be on the perfusion circle
    Segment s;
    s.myRadius = aRadius;
    s.myCoordinate = ptRoot;
    s.myIndex = 0;
    s.myKTerm = 0;
    myVectSegments.push_back(s);
    myVectParent.push_back(0); //if parent index is itsef it is the root (special segment of length 0).
    myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0)); // if children index is itself, it is an end segment.
    updateLengthFactor();
    
    // Construction of the first segment after the root
    Segment s1;
    s1.myRadius = aRadius;
    //s1.myCoordinate = generateRandomPtOnDisk(myTreeCenter, myRsupp);
    s1.myCoordinate = GeomHelpers::generateRandomPtOnDisk(myTreeCenter, my_rPerf);
    double myLength = (ptRoot-s1.myCoordinate).norm()*myLengthFactor;
    s1.myIndex = 1;
    s1.myKTerm = 1; //it contains terminal itself
    s1.myResistance = 8.0*my_mu*myLength/M_PI;
    s1.myFlow = my_qTerm;
    s1.myBeta = 1.0;
    
    myVectSegments.push_back(s1);
    myVectTerminals.push_back(1);
    myVectParent.push_back(0); //if parent index is the root
    myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0)); // if children index is itself, it is an end segment.
    updateRootRadius();
    DGtal::trace.info() << "Construction initialized..." << std::endl;
  };
  */
  /**
   * @brief Constructor.
   * @brief It generates the first root segment from a given terminal point.
   * @param ptCenter: coordinates of the center of the tree
   * @param ptRoot: coordinates of the root special segment (who have no parent)
   * @param ptTerm: coordinates of the first terminal point
   * @param nTerm: number of terminal segments.
   **/
  /*
  CoronaryArteryTree(const TPointD &ptCenter, const TPointD &ptRoot, const TPointD &ptTerm, unsigned int nTerm, int nNeighbor = 20, double aRadius = 1.0 ){
    assert(nTerm>=1);
    myTreeCenter = ptCenter;
    my_rPerf = (ptCenter - ptRoot).norm();
    if(TDim==2) {
      my_aPerf = M_PI*my_rPerf*my_rPerf;
      myRsupp = sqrt(my_aPerf/(nTerm*M_PI));
    }
    else {//TDim==3
      my_aPerf = 4.0*M_PI*my_rPerf*my_rPerf*my_rPerf/3.0;
      myRsupp = pow(3.0*my_aPerf/(4.0*M_PI*nTerm),1.0/3.0);
    }
    my_NTerm = nTerm;
    my_qTerm = my_qPerf / my_NTerm;
    if(nTerm > 250)
      my_pTerm = 7.98e3;
    my_pDrop = my_pPerf-my_pTerm;
    
    //myDThresold = sqrt(M_PI*myRsupp*myRsupp/myKTerm);
    
    // Construction of the special root segment
    assert((ptRoot - myTreeCenter).norm() == my_rPerf); //ptRoot must be on the perfusion circle
    Segment s;
    s.myRadius = aRadius;
    s.myCoordinate = ptRoot;
    s.myIndex = 0;
    s.myKTerm = 0;
    s.myFlow = my_qTerm;
    myVectSegments.push_back(s);
    myVectParent.push_back(0); //if parent index is itsef it is the root (special segment of length 0).
    myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0)); // if children index is itself, it is an end segment.
    updateLengthFactor();
    
    // Construction of the first segment after the root
    assert((ptTerm - myTreeCenter).norm() <= my_rPerf); //ptTerm must be in the perfusion disk
    Segment s1;
    s1.myRadius = aRadius;
    s1.myCoordinate = ptTerm;
    double myLength = (ptRoot-s1.myCoordinate).norm()*myLengthFactor;
    s1.myIndex = 1;
    s1.myKTerm = 1; //it contains terminal itself
    s1.myResistance = 8.0*my_mu*myLength/M_PI;
    s1.myFlow = my_qTerm;
    s1.myBeta = 1.0;
    
    myVectSegments.push_back(s1);
    myVectTerminals.push_back(1);
    myVectParent.push_back(0); //if parent index is the root
    myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0)); // if children index is itself, it is an end segment.
    updateRootRadius();
    DGtal::trace.info() << "Construction initialized..." << std::endl;
  };
  */
    /**
     * Constructor for base tests (used on various files of directory tests, like testGeom).
     * @param r: the radius of the domain circle
     **/
    
    CoronaryArteryTree(double r=1.0){
      myTreeCenter = TPointD::diagonal(0);
      myRsupp = r;
      my_aPerf = 1.0;
      my_NTerm = 1;
      // Construction of the special root segment
      Segment s;
      s.myRadius = 1.0;
      s.myCoordinate = TPointD(0,r);//ptRoot;
      s.myIndex = 0;
      myVectSegments.push_back(s);
      myVectParent.push_back(0); //if parent index is itsef it is the root (special segment of length 0).
      myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0)); // if children index is itself, it is an end segment.
      my_rPerf = myRsupp;
      
      // Construction of the first segment after the root
      Segment s1;
      s1.myRadius = 1.0;
      s1.myCoordinate = TPointD::diagonal(0);
      s1.myIndex = 1;
      myVectSegments.push_back(s1);
      myVectTerminals.push_back(1);
      myVectParent.push_back(0); //if parent index is the root
      myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0)); // if children index is itself, it is an end segment.
      DGtal::trace.info() << "Construction initialized..." << std::endl;
    };
    
    
  // ----------------------- Interface --------------------------------------
    
  /**
   * Tries to add a new segment from a given point and a nearest segment given by index.
   * @param p the extremity of the new segment to be created
   * @param segIndex the index of the near segment to p.
   * @param nbIter maximal number of iteration
   * @param tolerance convergence for tree volume gradient
   * @param verbose used to display process algorithm  informations.
   * @return true of the new segment is created, false in the other case.
   * (for instance if an intersection to previous point was present)
   **/
  bool isAddable(const TPointD &p,
                 unsigned int segIndex,
                 unsigned int nbIter,
                 double tolerance,
                 bool verbose = true);
  
  /**
   * Verifies if there is an intersection between a thick segment defined by two points
   * and all other segments of the tree
   * @param ptA the first point
   * @param ptA the second point
   * @param r : radius of the thick segment
   */
  bool isIntersectingTree(const TPointD &ptA,
                          const TPointD &ptB,
                          double r,
                          unsigned int idSeg) const;
  
  /**
   * Verifies if there is an intersection between a thick segment defined by two points
   * and all other segments of the tree, except the segments given in idExcept
   * @param ptA the first point
   * @param ptA the second point
   * @param r : radius of the thick segment
   * @param idExcept: indices of three segments to ignore
   */
  bool isIntersectingTree(const TPointD &ptA,
                          const TPointD &ptB,
                          double r,
                          std::tuple<int, int, int> idExcept) const;
  
  /**
   * Update the distribution of segmental flows after adding a new segment (new bifurcation)
   * @param segIndex index of the parent segment to be updated
   */
  CoronaryArteryTree::Segment updateResistanceFromRoot(unsigned int segIndex=1);
  
  /**
   * Updates the flow of a segment (as given in Eq. 11)
   * @param segIndex index of the parent segment to be updated
   */
  void updateFlow(unsigned int segIndex);
  
  /**
   * Updates the length factor value (as given in Eq. 9)
   */
  void updateLengthFactor();

  /**
   * Computes the distance threshold  (as given in Eq. 12)
   * @return the threshold
   */
  double getDistanceThreshold();
  
  /**
   * Updates the radius of a segment (as given in Eq. 19)
   * @param segIndex index of the parent segment to be updated
   * @param beta beta value of the segment
   */
  void updateRadius(unsigned int segIndex, double beta);
  
  /**
   * Updates the root radius after updating flow parameters (as given in Eq. 18)
   */
  void updateRootRadius();
  
  /**
   * Computes the total tree volume (as given in Eq. 20)
   */
  double computeTotalVolume(unsigned int segIndex = 1);
  
  /**
   * Computes the nearest segment of the given point
   * @param pt: a point
   * @return the index of the nearest segement
   */
  unsigned int getNearestSegment(const TPointD &pt);

  /**
   * Computes the barycenter of the given point and the two points of the segment
   * @param p: a point
   * @param segIndex : the index of the segement used for the computation
   */
  TPointD findBarycenter(const TPointD &p, unsigned int index);

  /**
   * From a segment returns a vector of segment index representing the path to the root.
   */
  std::vector<unsigned int> getPathToRoot(const Segment &s);
  
  /**
   * Generates a new location with distance constraints.
   * @param nbTrials : number of trials before reducing the distance constaint value
   *
   */
  
  DGtal::PointVector<TDim, double>
  generateNewLocation(unsigned int nbTrials = 1000);

  /**
   * Generates a new location with distance constraints.
   * @param myDThresold : the distance constaint value
   * @return the generated point and a bool true if sucess and false othewise
   */

  std::pair<DGtal::PointVector<TDim, double>, bool> generateALocation(double myDThresold);

  /**
   * Computes the length of a segment represented with the index and mutiliplies by the length factor.
   * @param segIndex : the index of the segement used for the comparison
   */
  double getLengthSegment(unsigned int segIndex);
  
  /**
   * Computes the projected distance from a segment represented with the index  and the point given as argument.
   * @param index : the index of the segement used for the comparison
   * @param p : a point
   */
  double getProjDistance(unsigned int index, const TPointD &p) const;
  /**
   * Computes the projected distance from a segment represented with the index  and the point given as argument.
   * @param p0 : a point representing one extremity
   * @param p1 : a point representing another extremity
   * @param p : a point to be projected
   */
  double getProjDistance(const TPointD &p0, const TPointD &p1, const TPointD &p) const;
  
  /**
   * Check if a new added point is too close the nearest segment.
   * @param p : a point
   * @param minDist: min distance
   * @return true is a point is too close
   */
  bool isToCloseFromNearest(const TPointD &p, double minDist) const;
  
  /**
   * Computes the n nearest index of neighborhood segments  of a given point
   * @param p : the point considered to get the nearest point
   * @param n : the number of nearest point to be recovered
   */
  std::vector<unsigned int> getN_NearestSegments(const TPointD &p,
                                                 unsigned int n) const;
  
  /**
   * Verifies the condition of degenerate after Kamyia solution : 2 ri < li
   * @param length length segment
   * @param radius radius segment
   * @return true if degenerate and false otherwise
   */
  bool isDegenerate(double length, double radius);
  
  /**
   * Solver for  Kamyia optimization
   * @param pCurrent: the coordinate of the current bifurcation position
   * @param pParent: the coordinate of the parent point (the point from which the bifucation arrives)
   * @param sCurrent: the midle segment (from pParent point)
   * @param sL: the left segment
   * @param sR: the right segment
   * @param pOpt: output opmized point
   * @param r0: output radius of the middle segment
   * @param r1: output radius of the left segment
   * @param r2: output radius of the right segment
   */
  bool kamyiaOptimization(const TPointD& pCurrent,
                          const TPointD& pParent,
                          double rCurrent,
                          const Segment &sL,
                          const Segment &sR,
                          unsigned int nbIter,
                          TPointD& pOpt,
                          double& r0,
                          double& r1,
                          double& r2);
  
  /**
   * Use to display object.
   * @param out  the object of class 'GeodesicGraphComputer' to write.
   */
  void selfDisplay( std::ostream & out ) const;
  
  /**
   * Get the supported domain of the tree. By default it is defined from the circle center.
   * If the domain if defined from a mask image, the center if computed from the imate center.
   */
  TPoint getDomainCenter() const;
  
  
  /**
   * Export the current display of the tree.
   */
  
  void boardDisplay(double thickness = 1,
                    bool clearDisplay = true);
  void exportBoardDisplay(const std::string &fileName = "result.eps",
                          double thickness = 1,
                          bool updateDisplay = true,
                          bool clearDisplay = true);
  
  /**
   * Restrain the domain from the image mask.
   * @param image name that should represent the restricted domain (ie all pixels set to 255)
   * @param threshold the threshold to consider foreground value (used to test if at least one domain pixel exist).
   * @return true the restriction was well set (image exist).
   */
  bool restrainDomain(const Image &imageDom, unsigned int threshold=128);

  
  /**
   * Try to find the root point (the non distal vertex of the first segment).
   * The center of tree (distal segment) is supposed to be already given and
   * and this function try to find a point P_root on circular area with the condition that
   * that the segment [P_root, P_center] .
   * @param d the distance to search around the central point.
   *  @param[out] ptRoot the root point is updated if found.
   * @return true of the root point was found.
   */
   bool searchRootFarthest(const double & d, TPointD &ptRoot );
  
  
private:

  
  
};



/**
 * Overloads 'operator<<' for displaying objects of class 'CoronaryArteryTree'.
 * @param out the output stream where the object is written.
 * @param aCoronaryTree the object of class 'CoronaryArteryTree' to write.
 * @return the output stream after the writing.
 */
template <int TDim>
std::ostream&
operator<< ( std::ostream & out, const CoronaryArteryTree<TDim> & aCoronaryTree );

template <int TDim>
bool
operator==(typename CoronaryArteryTree<TDim>::Segment  S1,typename  CoronaryArteryTree<TDim>::Segment S2);

template <int TDim>
bool
operator!=(typename CoronaryArteryTree<TDim>::Segment  S1, typename CoronaryArteryTree<TDim>::Segment S2);







///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods 
#include "CoronaryArteryTree.ih"



///////////////////////////////////////////////////////////////////////////////

#endif // !defined CORONARY_ARTERY_TREE_h

#undef CORONARY_ARTERY_TREE_RECURSES
#endif // else defined(CORONARY_ARTERY_TREE_RECURSES)





