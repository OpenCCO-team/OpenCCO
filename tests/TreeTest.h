#pragma once

#if defined(TREE_TEST_RECURSES)
#error Recursive header files inclusion detected in CoronaryArteryTree.h
#else // defined(TREE_TEST_RECURSES)
/** Prevents recursive inclusion of headers. */
#define TREE_TEST_RECURSES

#if !defined TREETEST_h
/** Prevents repeated inclusion of headers. */
#define TREETEST_h



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








class TreeTest: public CoronaryArteryTree<ImageMaskDomainCtrl<2>, 2>
{
public:
    // Constructor do nothing mainly used for specific test
    TreeTest(ImageMaskDomainCtrl<2> &ctr, double r=1.0 ): CoronaryArteryTree(ctr){
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
    }
    /*
     * Constructor used mainly in testCompCCO
     */
    
    TreeTest(const TPointD &ptCenter, const TPointD &ptRoot, const TPointD &ptTerm, unsigned int nTerm,  ImageMaskDomainCtrl<2> &ctr , double aRadius = 1.0 ): CoronaryArteryTree(ctr){
        assert(nTerm>=1);
        myTreeCenter = ptCenter;
        my_rPerf = (ptCenter - ptRoot).norm();
        my_aPerf = M_PI*my_rPerf*my_rPerf;
        myRsupp = sqrt(my_aPerf/(my_NTerm*M_PI));
        my_NTerm = my_NTerm;
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


   /**
   * @brief Constructor.
   * @brief It generates the first root segment by randomly choose the first terminal point.
   * @param ptRoot: coordinates of the root special segment (who have no parent)
   * @param aPerf: surface of the perfusion.
   * @param nTerm: number of terminal segments.
   **/
  
  TreeTest(const TPointD &ptRoot, double aPerf, unsigned int nTerm,
           ImageMaskDomainCtrl<2> &ctr,double aRadius = 1.0 ): CoronaryArteryTree(ctr){
    assert(nTerm>=1);
    myTreeCenter = TPointD::diagonal(0.0);
    myRsupp = sqrt(aPerf/(nTerm*M_PI));
    my_rPerf = sqrt(aPerf/M_PI);
    
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
    s1.myCoordinate = myDomainController().randomPoint();
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
  
    bool
    addSegmentFromPoint(const TPointD &p,
                        unsigned int nearIndex)
    {
      TPointD newCenter = findBarycenter(p, nearIndex);
      bool inter = hasNearestIntersections(myVectParent[nearIndex], nearIndex, p, newCenter,  10);
      double minDistance = 5.0;
      if (inter){
        DGtal::trace.warning() << "detection intersection" << std::endl;
        return false;
      }
      // Check barycenter is not too close considered segment
      if ((newCenter-myVectSegments[nearIndex].myCoordinate).norm() < minDistance ||
          (newCenter-myVectSegments[myVectParent[nearIndex]].myCoordinate).norm() < minDistance ||
          (newCenter - p).norm() < minDistance){
        DGtal::trace.warning() << "new barycenter too close!!!!!!!!" << std::endl;
        return false;
      }
      // Check new point with new segment of barycenter:
      // - newPt and new segment [Barycenter-OriginNewSeg]
      // - newPt and new segment [Barycenter-FatherNewSeg]
      if ((newCenter-myVectSegments[nearIndex].myCoordinate).norm() < minDistance ||
          (newCenter-myVectSegments[myVectParent[nearIndex]].myCoordinate).norm() < minDistance) {
        DGtal::trace.warning() << "initial too close to new!!!!!!!!" << std::endl;
        return false;
      }
      if (getProjDistance(nearIndex, p) < minDistance) {
        DGtal::trace.warning() << "initial too close!!!!!!!!" << std::endl;
        return false;

      }
      if (isToCloseFromNearest(p, minDistance)){
        DGtal::trace.warning() << "detection near too close " << std::endl;
        return false;
      }
      if (getProjDistance(nearIndex, p) < minDistance||
          getProjDistance(myVectParent[nearIndex], p) < minDistance||
          getProjDistance(myVectChildren[nearIndex].first, p) < minDistance||
          getProjDistance(myVectChildren[nearIndex].second, p) < minDistance){
        DGtal::trace.warning() << "detection too close existing" << std::endl;
        return false;
      }
      // Creation of the left child
      Segment sNewLeft;
      sNewLeft.myCoordinate = myVectSegments[nearIndex].myCoordinate;
      sNewLeft.myIndex = (unsigned int) myVectSegments.size();
      sNewLeft.myRadius = myVectSegments[nearIndex].myRadius;
      sNewLeft.myFlow = myVectSegments[nearIndex].myFlow;
      sNewLeft.myKTerm = myVectSegments[nearIndex].myKTerm;
      sNewLeft.myResistance = myVectSegments[nearIndex].myResistance;
      myVectSegments.push_back(sNewLeft);
      myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0));
      myVectParent.push_back(nearIndex);
      
      // Update for the new parent of the new left segment
      unsigned int leftGrandChildIndex = myVectChildren[myVectSegments[nearIndex].myIndex].first;
      unsigned int rightGrandChildIndex = myVectChildren[myVectSegments[nearIndex].myIndex].second;
      myVectParent[leftGrandChildIndex] = sNewLeft.myIndex;
      myVectParent[rightGrandChildIndex] = sNewLeft.myIndex;

      // Update of the child of the new left segment (sNewLeft)
      myVectChildren[sNewLeft.myIndex].first = leftGrandChildIndex;
      myVectChildren[sNewLeft.myIndex].second = rightGrandChildIndex;

      // Creation of the right child
      Segment sNewRight;
      sNewRight.myCoordinate = p;
      sNewRight.myRadius = myVectSegments[nearIndex].myRadius;
      sNewRight.myIndex = (unsigned int) myVectSegments.size();
      sNewRight.myFlow = my_qTerm;
      sNewRight.myKTerm = 1;
      myVectSegments.push_back(sNewRight);
      myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0));
      myVectParent.push_back(nearIndex);
      myVectTerminals.push_back(sNewRight.myIndex);
      
      // Update center segment
      myVectSegments[nearIndex].myCoordinate = newCenter;
      myVectSegments[nearIndex].myFlow = myVectSegments[nearIndex].myFlow + my_qTerm;
      myVectSegments[nearIndex].myKTerm = sNewLeft.myKTerm + 1;
      //update childrens of center segment
      myVectChildren[nearIndex].first = sNewLeft.myIndex;
      myVectChildren[nearIndex].second = sNewRight.myIndex;
      
      // Update physilogique paramaters
        myKTerm++;
      updateFlow(myVectParent[nearIndex]);
      updateResistanceFromRoot();
      updateRootRadius();
      
      return true;
    }

    
    
};



class TreeTestCirc: public CoronaryArteryTree<CircularDomainCtrl<2>, 2>
{
public:
    // Constructor do nothing mainly used for specific test
    TreeTestCirc(CircularDomainCtrl<2> &ctr, double r=1.0): CoronaryArteryTree(ctr){
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
    }
    /*
     * Constructor used mainly in testCompCCO
     */
    
    TreeTestCirc(const TPointD &ptCenter, const TPointD &ptRoot, const TPointD &ptTerm, unsigned int nTerm,  CircularDomainCtrl<2> &ctr, double aRadius = 1.0 ): CoronaryArteryTree(ctr){
        assert(nTerm>=1);
        myTreeCenter = ptCenter;
        my_rPerf = (ptCenter - ptRoot).norm();
        my_aPerf = M_PI*my_rPerf*my_rPerf;
        myRsupp = sqrt(my_aPerf/(my_NTerm*M_PI));
        my_NTerm = my_NTerm;
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
    
    TreeTestCirc(const TPointD &ptRoot, double aPerf, unsigned int nTerm,
                 CircularDomainCtrl<2> &ctr, double aRadius = 1.0 ): CoronaryArteryTree(ctr){
      assert(nTerm>=1);
      myTreeCenter = TPointD::diagonal(0.0);
      myRsupp = sqrt(aPerf/(nTerm*M_PI));
      my_rPerf = sqrt(aPerf/M_PI);
      
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
      s1.myCoordinate = myDomainController().randomPoint();
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
    bool
    addSegmentFromPoint(const TPointD &p,
                        unsigned int nearIndex)
    {
      TPointD newCenter = findBarycenter(p, nearIndex);
      bool inter = hasNearestIntersections(myVectParent[nearIndex], nearIndex, p, newCenter,  10);
      double minDistance = 5.0;
      if (inter){
        DGtal::trace.warning() << "detection intersection" << std::endl;
        return false;
      }
      // Check barycenter is not too close considered segment
      if ((newCenter-myVectSegments[nearIndex].myCoordinate).norm() < minDistance ||
          (newCenter-myVectSegments[myVectParent[nearIndex]].myCoordinate).norm() < minDistance ||
          (newCenter - p).norm() < minDistance){
        DGtal::trace.warning() << "new barycenter too close!!!!!!!!" << std::endl;
        return false;
      }
      // Check new point with new segment of barycenter:
      // - newPt and new segment [Barycenter-OriginNewSeg]
      // - newPt and new segment [Barycenter-FatherNewSeg]
      if ((newCenter-myVectSegments[nearIndex].myCoordinate).norm() < minDistance ||
          (newCenter-myVectSegments[myVectParent[nearIndex]].myCoordinate).norm() < minDistance) {
        DGtal::trace.warning() << "initial too close to new!!!!!!!!" << std::endl;
        return false;
      }
      if (getProjDistance(nearIndex, p) < minDistance) {
        DGtal::trace.warning() << "initial too close!!!!!!!!" << std::endl;
        return false;

      }
      if (isToCloseFromNearest(p, minDistance)){
        DGtal::trace.warning() << "detection near too close " << std::endl;
        return false;
      }
      if (getProjDistance(nearIndex, p) < minDistance||
          getProjDistance(myVectParent[nearIndex], p) < minDistance||
          getProjDistance(myVectChildren[nearIndex].first, p) < minDistance||
          getProjDistance(myVectChildren[nearIndex].second, p) < minDistance){
        DGtal::trace.warning() << "detection too close existing" << std::endl;
        return false;
      }
      // Creation of the left child
      Segment sNewLeft;
      sNewLeft.myCoordinate = myVectSegments[nearIndex].myCoordinate;
      sNewLeft.myIndex = (unsigned int) myVectSegments.size();
      sNewLeft.myRadius = myVectSegments[nearIndex].myRadius;
      sNewLeft.myFlow = myVectSegments[nearIndex].myFlow;
      sNewLeft.myKTerm = myVectSegments[nearIndex].myKTerm;
      sNewLeft.myResistance = myVectSegments[nearIndex].myResistance;
      myVectSegments.push_back(sNewLeft);
      myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0));
      myVectParent.push_back(nearIndex);
      
      // Update for the new parent of the new left segment
      unsigned int leftGrandChildIndex = myVectChildren[myVectSegments[nearIndex].myIndex].first;
      unsigned int rightGrandChildIndex = myVectChildren[myVectSegments[nearIndex].myIndex].second;
      myVectParent[leftGrandChildIndex] = sNewLeft.myIndex;
      myVectParent[rightGrandChildIndex] = sNewLeft.myIndex;

      // Update of the child of the new left segment (sNewLeft)
      myVectChildren[sNewLeft.myIndex].first = leftGrandChildIndex;
      myVectChildren[sNewLeft.myIndex].second = rightGrandChildIndex;

      // Creation of the right child
      Segment sNewRight;
      sNewRight.myCoordinate = p;
      sNewRight.myRadius = myVectSegments[nearIndex].myRadius;
      sNewRight.myIndex = (unsigned int) myVectSegments.size();
      sNewRight.myFlow = my_qTerm;
      sNewRight.myKTerm = 1;
      myVectSegments.push_back(sNewRight);
      myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0));
      myVectParent.push_back(nearIndex);
      myVectTerminals.push_back(sNewRight.myIndex);
      
      // Update center segment
      myVectSegments[nearIndex].myCoordinate = newCenter;
      myVectSegments[nearIndex].myFlow = myVectSegments[nearIndex].myFlow + my_qTerm;
      myVectSegments[nearIndex].myKTerm = sNewLeft.myKTerm + 1;
      //update childrens of center segment
      myVectChildren[nearIndex].first = sNewLeft.myIndex;
      myVectChildren[nearIndex].second = sNewRight.myIndex;
      
      // Update physilogique paramaters
        myKTerm++;
      updateFlow(myVectParent[nearIndex]);
      updateResistanceFromRoot();
      updateRootRadius();
      
      return true;
    }

    
    
};


///////////////////////////////////////////////////////////////////////////////

#endif // !defined TREETEST_h

#undef TREE_TEST_RECURSES
#endif // else defined(TREE_TEST_RECURSES)


