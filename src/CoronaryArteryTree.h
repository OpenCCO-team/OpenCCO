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

#if defined(CORONARY_ARTERY_TREE_RECURSES)
#error Recursive header files inclusion detected in CoronaryArteryTree.h
#else // defined(CORONARY_ARTERY_TREE_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CORONARY_ARTERY_TREE_RECURSES

#if !defined CORONARY_ARTERY_TREE_h
/** Prevents repeated inclusion of headers. */
#define CORONARY_ARTERY_TREE_h

#include <type_traits>


#include <iostream>
#include <math.h>
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/images/ImageSelector.h"

#include "DGtal/helpers/StdDefs.h"

#include "GeomHelpers.h"
#include "DomainController.h"

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
 * @brief Class to represent and construct a Coronary Artery Tree used through the CCO algorithm.
 *
 *
 */

template <class DomCtr, int TDim>
class CoronaryArteryTree{
    /**
     * Class to handle the representation of the coronary tree.
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
        // hydrodynamc resistance (R star)
        double myResistance = 0.00;
        // flow (Qi)
        double myFlow = 0.00;
        // radius ratio of the segment and its parent
        double myBeta = 1.0;
    };
    
    // to recover the Children left (first) and right (second) on an indexed segment.
    std::vector< SegmentChildren >  myVectChildren;
    // to recover the parent of an indexed segment
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
        
        // my_pTerm pressure at distal ends of terminal segment
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
        unsigned int myNumNeighbor = 20;
        
        // myLengthFactor : scale factor (updated during tree growth after each added bifurcation)
        double myLengthFactor = 1.0;
    // End: Internal algorithm parameter
    //-----------------------------
    
    
   
    std::reference_wrapper<DomCtr> myDomainController_;
    
    DomCtr& myDomainController() const {
        return myDomainController_.get();
    };
    
    
    // To handle display
    DGtal::Board2D myBoard;
    
    
    
public:
    
    // Constructor do nothing mainly used for specific test
    CoronaryArteryTree(DomCtr &aDomCtr): myDomainController_(aDomCtr){}
    /**
     * @brief Default constructor.
     * @brief It generates the first root segment with tree center as the first terminal point.
     * @param aPerf: surface of the perfusion.
     * @param nTerm: number of terminal segments.
     **/
    
    CoronaryArteryTree(double aPerf, unsigned int nTerm,  DomCtr &aDomCtr,
                       double aRadius = 1.0): myDomainController_(aDomCtr) {
        assert(nTerm>=1);
        
        for (auto i=0; i < TDim; i++){myTreeCenter[i]=aDomCtr.myCenter[i];}
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
        myVectParent.push_back(0); // parent index is the root
        myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0)); // if children index is itself, it is an end segment.
        updateRootRadius();
        
        // special case for implicit domain:
        if (aDomCtr.myUpdateType == DomCtr::UPDATED )
        {
            aDomCtr.myRadius = my_rPerf;
        }
        if(aDomCtr.randomPoint() == TPoint::diagonal(0) )
        {
            DGtal::trace.error() << "Domain too restrained, not possible to find random"
            << " candidate in domain (probably reduce the minimal distance to border)"
            << std::endl;
            exit(1);
        }

        DGtal::trace.info() << "Construction initialized..." << std::endl;

        myBoard.setLineCap(LibBoard::Shape::LineCap::RoundCap);
    };
    
    
   
    
    
    // ----------------------- Interface --------------------------------------
    
  
    /**
     * Tries to add a new segment from a given point and a nearest segment given by index.
     * @param p the extremity of the new segment to be created
     * @param segIndex the index of the near segment to p.
     * @param nbIter maximal number of iteration
     * @param tolerance convergence boundary for tree volume gradient
     * @param nbNeibour number of neighbours to be considered for intersecting test
     * @param verbose used to display process algorithm  informations.
     * @return true if the new segment is created, false in the other case.
     * (for instance if an intersection to previous point was present)
     **/
    bool isAddable(const TPointD &p,
                   unsigned int segIndex,
                   unsigned int nbIter,
                   double tolerance,
                   unsigned int nbNeibour = 10,
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
     * @return the index of the nearest segment
     */
    unsigned int getNearestSegment(const TPointD &pt);
    
    /**
     * Computes the barycenter of the given point and the two points of the segment
     * @param p: a point
     * @param segIndex : the index of the segment used for the computation
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
     * @param segIndex : the index of the segment used for the comparison
     */
    double getLengthSegment(unsigned int segIndex);
    
    /**
     * Computes the projected distance from a segment represented with the index  and the point given as argument.
     * @param index : the index of the segment used for the comparison
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


    void createPointInDomain(int nbpoints);

    void exportBoardDisplay(const std::string &fileName = "result.eps",
                            double thickness = 1,
                            bool updateDisplay = true,
                            bool clearDisplay = true);
    
    
    /**
     * Try to find the root point (the non distal vertex of the first segment).
     * The center of tree (distal segment) is supposed to be already given and
     * and this function try to find a point P_root on circular area with the condition that
     * that the segment [P_root, P_center] .
     * @param d the distance to search around the central point.
     *  @param[out] ptRoot the root point is updated if found.
     * @return true of the root point was found.
     */
    //   bool searchRootFarthest(const double & d, TPointD &ptRoot );
    //   // 3d specialisation
    //   bool searchRootFarthest3d(const double & d, TPointD &ptRoot );
    //   // 2d specialisation
    //   bool searchRootFarthest2d(const double & d, TPointD &ptRoot );
    
    
private:
    
    
    
};


/**
 * Overloads 'operator<<' for displaying objects of class 'CoronaryArteryTree'.
 * @param  out the output stream where the object is written.
 * @param aCoronaryTree the object of class 'CoronaryArteryTree' to write.
 * @return the output stream after the writing.
 */
template <typename DomCtr, int TDim>
std::ostream&
operator<< ( std::ostream & out, const CoronaryArteryTree<DomCtr, TDim> & aCoronaryTree );

template <typename DomCtr, int TDim>
bool
operator==(typename CoronaryArteryTree<DomCtr, TDim>::Segment  S1,typename  CoronaryArteryTree<DomCtr, TDim>::Segment S2);

template <typename DomCtr, int TDim>
bool
operator!=(typename CoronaryArteryTree<DomCtr, TDim>::Segment  S1, typename CoronaryArteryTree<DomCtr, TDim>::Segment S2);







///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods 
#include "CoronaryArteryTree.ih"



///////////////////////////////////////////////////////////////////////////////

#endif // !defined CORONARY_ARTERY_TREE_h

#undef CORONARY_ARTERY_TREE_RECURSES
#endif // else defined(CORONARY_ARTERY_TREE_RECURSES)


