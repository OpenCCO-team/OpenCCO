
#include "CoronaryArteryTree.h"
#include "DGtal/io/boards/Board2D.h"
#include "geomhelpers.h"
#include "DGtal/io/readers/PPMReader.h"
#include "DGtal/images/ArrayImageIterator.h"
#include "DGtal/io/readers/PGMReader.h"
#include <math.h>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/io/colormaps/GradientColorMap.h"


using namespace std;
bool
CoronaryArteryTree::addSegmentFromPoint(const Point2D &p)
{
  unsigned int nearIndex = getNearestSegment(p);
  return addSegmentFromPoint(p, nearIndex);
}

void
CoronaryArteryTree::updateResistanceTerminal(unsigned int segIndex)
{
  assert(myVectChildren[segIndex].first==0);
  assert(myVectChildren[segIndex].second==0);
  //Update HydroResistance
  double l = getLengthSegment(segIndex);
  myVectSegments[segIndex].myResistance = 8.0*my_nu*l/M_PI;
  //myVectSegments[segIndex].myResistance = 8.0*my_nu*getLengthSegment(segIndex)/(M_PI*pow(myVectSegments[segIndex].myRadius,4));
  assert(myVectSegments[segIndex].myResistance != 0);
  //std::cout<<"myResistance="<<myVectSegments[segIndex].myResistance<<std::endl;
}

void
CoronaryArteryTree::updateResistance(unsigned int segIndex)
{
  if(myVectSegments[segIndex].myIndex != 0) {
    //Update HydroResistance
    double l = getLengthSegment(segIndex);
    myVectSegments[segIndex].myResistance = 8.0*my_nu*l/M_PI;
    //myVectSegments[segIndex].myResistance = 8.0*my_nu*getLengthSegment(segIndex)/(M_PI*pow(myVectSegments[segIndex].myRadius,4));
    if(myVectChildren[segIndex].first!=0 && myVectChildren[segIndex].second!=0) {
      // Get left and right children of the segment
      Segment<Point2D> sLeft = myVectSegments[myVectChildren[segIndex].first];
      Segment<Point2D> sRight = myVectSegments[myVectChildren[segIndex].second];
      double r1 = sLeft.myBeta; //(sLeft.myRadius/myVectSegments[segIndex].myRadius);
      double r2 = sRight.myBeta; //(sRight.myRadius/myVectSegments[segIndex].myRadius);
      //double rr1 = (r1*r1*r1*r1)/sLeft.myResistance;
      //double rr2 = (r2*r2*r2*r2)/sRight.myResistance;
      myVectSegments[segIndex].myResistance += 1.0/((r1*r1*r1*r1)/sLeft.myResistance + (r2*r2*r2*r2)/sRight.myResistance) ;
    }
    assert(myVectSegments[segIndex].myResistance != 0);
    //std::cout<<"myResistance="<<myVectSegments[segIndex].myResistance<<std::endl;
    
    std::pair<int, int> children = myVectChildren[myVectParent[segIndex]];
    int brotherIndex;
    if(children.first != segIndex)
      brotherIndex = children.first;
    else
      brotherIndex = children.second;
    //Update brother resistance
    if(brotherIndex<segIndex)
      updateResistance(brotherIndex);
    else {
      //Update beta
      if(myVectSegments[myVectParent[segIndex]].myIndex!=0) { //it's not the first segment
        double segFlow = myVectSegments[segIndex].myFlow;
        double brotherFlow = myVectSegments[brotherIndex].myFlow;
        double segHydro = myVectSegments[segIndex].myResistance;
        double brotherHydro = myVectSegments[brotherIndex].myResistance;
        double ratioR = pow((segFlow*segHydro) / (brotherFlow*brotherHydro), 0.25); //ratio of radii of two brother segments
        //ratio of radii of current segment wrt the parent segment
        myVectSegments[segIndex].myBeta = pow(1.0/(1+pow(ratioR,-my_gamma)),1.0/my_gamma);
        myVectSegments[brotherIndex].myBeta = pow(1.0/(1+pow(1.0/ratioR,-my_gamma)),1.0/my_gamma);
      }
    }
    
    //Update parent resistance
    updateResistance(myVectSegments[myVectParent[segIndex]].myIndex);
  }
}

void
CoronaryArteryTree::updateResistance(unsigned int segIndex, int order)
{
  if(myVectSegments[segIndex].myIndex != 0) {
    //Update HydroResistance
    double l = getLengthSegment(segIndex);
    myVectSegments[segIndex].myResistance = 8.0*my_nu*l/M_PI;
    //myVectSegments[segIndex].myResistance = 8.0*my_nu*getLengthSegment(segIndex)/(M_PI*pow(myVectSegments[segIndex].myRadius,4));
    if(myVectChildren[segIndex].first!=0 && myVectChildren[segIndex].second!=0) {
      // Get left and right children of the segment
      Segment<Point2D> sLeft = myVectSegments[myVectChildren[segIndex].first];
      Segment<Point2D> sRight = myVectSegments[myVectChildren[segIndex].second];
      double r1 = sLeft.myBeta; //(sLeft.myRadius/myVectSegments[segIndex].myRadius);
      double r2 = sRight.myBeta; //(sRight.myRadius/myVectSegments[segIndex].myRadius);
      //double rr1 = (r1*r1*r1*r1)/sLeft.myResistance;
      //double rr2 = (r2*r2*r2*r2)/sRight.myResistance;
      myVectSegments[segIndex].myResistance += 1.0/((r1*r1*r1*r1)/sLeft.myResistance + (r2*r2*r2*r2)/sRight.myResistance) ;
    }
    assert(myVectSegments[segIndex].myResistance != 0);
    //std::cout<<"myResistance="<<myVectSegments[segIndex].myResistance<<std::endl;
    
    std::pair<int, int> children = myVectChildren[myVectParent[segIndex]];
    int brotherIndex;
    if(children.first != segIndex)
      brotherIndex = children.first;
    else
      brotherIndex = children.second;
    //Update brother resistance
    if(order==0)
      updateResistance(brotherIndex, 1);
    else {
      //Update beta
      if(myVectSegments[myVectParent[segIndex]].myIndex!=0) { //it's not the first segment
        double segFlow = myVectSegments[segIndex].myFlow;
        double brotherFlow = myVectSegments[brotherIndex].myFlow;
        double segHydro = myVectSegments[segIndex].myResistance;
        double brotherHydro = myVectSegments[brotherIndex].myResistance;
        double ratioR = pow((segFlow*segHydro) / (brotherFlow*brotherHydro), 0.25); //ratio of radii of two brother segments
        //ratio of radii of current segment wrt the parent segment
        myVectSegments[segIndex].myBeta = pow(1.0/(1+pow(ratioR,-my_gamma)),1.0/my_gamma);
        myVectSegments[brotherIndex].myBeta = pow(1.0/(1+pow(1.0/ratioR,-my_gamma)),1.0/my_gamma);
      }
    }
    //Update parent resistance
    updateResistance(myVectSegments[myVectParent[segIndex]].myIndex, 0);
  }
}

CoronaryArteryTree::Segment<CoronaryArteryTree::Point2D>
CoronaryArteryTree::updateResistanceFromRoot(unsigned int segIndex) {
  if(segIndex != 0) {
    //Update HydroResistance
    double l = getLengthSegment(segIndex);
    myVectSegments[segIndex].myResistance = 8.0*my_nu*l/M_PI;
    if(myVectChildren[segIndex].first==0 && myVectChildren[segIndex].second==0)
      return myVectSegments[segIndex];
    
    std::pair<int, int> children = myVectChildren[segIndex];
    int segIndexLeft = children.first;
    int segIndexRight = children.second;
    //Update resitance of the two children segments
    Segment<Point2D> sLeft = updateResistanceFromRoot(segIndexLeft);
    Segment<Point2D> sRight = updateResistanceFromRoot(segIndexRight);
    
    //Update beta of the two children segments
    double segFlow = myVectSegments[segIndexLeft].myFlow;
    double brotherFlow = myVectSegments[segIndexRight].myFlow;
    double segHydro = myVectSegments[segIndexLeft].myResistance;
    double brotherHydro = myVectSegments[segIndexRight].myResistance;
    double ratioR = pow((segFlow*segHydro) / (brotherFlow*brotherHydro), 0.25); //ratio of radii of two brother segments
    //ratio of radii of current segment wrt the parent segment
    myVectSegments[segIndexLeft].myBeta = pow(1.0/(1+pow(ratioR,-my_gamma)),1.0/my_gamma);
    myVectSegments[segIndexRight].myBeta = pow(1.0/(1+pow(1.0/ratioR,-my_gamma)),1.0/my_gamma);
    //Compute resistance of parent segment
    double r1 = myVectSegments[segIndexLeft].myBeta;
    double r2 = myVectSegments[segIndexRight].myBeta;
    myVectSegments[segIndex].myResistance += 1.0/((r1*r1*r1*r1)/sLeft.myResistance + (r2*r2*r2*r2)/sRight.myResistance) ;
    return myVectSegments[segIndex];
  }
}
/*
void
CoronaryArteryTree::updateBeta() {
  for(size_t i=1; i<myVectSegments.size(); i++) {
    //get children of the seg
    int segIndex = myVectSegments[i].myIndex;
    std::pair<int, int> children = myVectChildren[myVectParent[segIndex]];
    if(children.first !=0 && children.second!=0) { //if the segment has children
      int sLIndex = children.first;
      int sRIndex = children.second;
      double sLFlow = myVectSegments[sLIndex].myFlow;
      double sRFlow = myVectSegments[sRIndex].myFlow;
      double sLHydro = myVectSegments[sLIndex].myResistance;
      double sRHydro = myVectSegments[sRIndex].myResistance;
      double ratioR = pow((sLFlow*sLHydro) / (sRFlow*sRHydro), 0.25); //ratio of radii of two brother segments
      //ratio of radii of current segment wrt the parent segment
      myVectSegments[sLIndex].myBeta = pow(1.0/(1+pow(ratioR,-my_gamma)),1.0/my_gamma);
      myVectSegments[sRIndex].myBeta = pow(1.0/(1+pow(1.0/ratioR,-my_gamma)),1.0/my_gamma);
    }
  }
}
*/
/*
void
CoronaryArteryTree::updateBeta(unsigned int segIndex) {
  if(myVectSegments[myVectParent[segIndex]].myIndex!=0) { //it's not the first segment
    std::pair<int, int> children = myVectChildren[myVectParent[segIndex]];
    int brotherIndex;
    if(children.first != segIndex)
      brotherIndex = children.first;
    else
      brotherIndex = children.second;
    
    double segFlow = myVectSegments[segIndex].myFlow;
    double brotherFlow = myVectSegments[brotherIndex].myFlow;
    double segHydro = myVectSegments[segIndex].myResistance;
    double brotherHydro = myVectSegments[brotherIndex].myResistance;
    double ratioR = pow((segFlow*segHydro) / (brotherFlow*brotherHydro), 0.25); //ratio of radii of two brother segments
    //ratio of radii of current segment wrt the parent segment
    myVectSegments[segIndex].myBeta = pow(1.0/(1+pow(ratioR,-my_gamma)),1.0/my_gamma);
    myVectSegments[brotherIndex].myBeta = pow(1.0/(1+pow(1.0/ratioR,-my_gamma)),1.0/my_gamma);
    //updateBeta(myVectParent[segIndex]);
    //updateBeta();
  }
}
*/
/*
void
CoronaryArteryTree::updateKTerm(unsigned int segIndex)
{
  Segment<Point2D> sC = myVectSegments[segIndex];
  while (sC.myIndex != 0){
    myVectSegments[segIndex].myKTerm ++;
    sC = myVectSegments[myVectParent[sC.myIndex]];
  }
}
*/
void
CoronaryArteryTree::updateFlow(unsigned int segIndex)
{
  Segment<Point2D> sC = myVectSegments[segIndex];
  while (sC.myIndex != 0){
    myVectSegments[sC.myIndex].myKTerm ++;
    //Eq. 11
    myVectSegments[sC.myIndex].myFlow = myVectSegments[sC.myIndex].myKTerm*my_qTerm;
    sC = myVectSegments[myVectParent[sC.myIndex]];
  }
}
/*
void
CoronaryArteryTree::updateFlow()
{
  for(size_t it=0; it<myVectSegments.size(); it++) {
    if(myVectChildren[it].first==0 && myVectChildren[it].second==0)
      myVectSegments[it].myFlow = my_qTerm;
    else //Eq. 11
      myVectSegments[it].myFlow = myVectSegments[it].myKTerm*my_qTerm;
  }
}
*/
void
CoronaryArteryTree::updateLengthFactor()
{
  double r_pk = sqrt(myKTerm + 1)*myRsupp; //Eq 9 => but not the same
  myLengthFactor =  r_pk / my_rPerf;
}

void
CoronaryArteryTree::updateRadius(unsigned int segIndex)
{
  if(myVectChildren[segIndex].first!=0 && myVectChildren[segIndex].second!=0) {
    //update radii left and right children
    Segment<Point2D> sLeft = myVectSegments[myVectChildren[segIndex].first];
    Segment<Point2D> sRight = myVectSegments[myVectChildren[segIndex].second];
    //Eq 19
    double rL = myVectSegments[segIndex].myRadius * sLeft.myBeta;
    myVectSegments[myVectChildren[segIndex].first].myRadius = rL;
    double rR = myVectSegments[segIndex].myRadius * sRight.myBeta;
    myVectSegments[myVectChildren[segIndex].second].myRadius = rR;
    updateRadius(sLeft.myIndex);
    updateRadius(sRight.myIndex);
  }
}

void
CoronaryArteryTree::updateRadius(unsigned int segIndex, double beta)
{
  if(myVectChildren[segIndex].first!=0 && myVectChildren[segIndex].second!=0) {
    //update radii left and right children
    Segment<Point2D> sLeft = myVectSegments[myVectChildren[segIndex].first];
    Segment<Point2D> sRight = myVectSegments[myVectChildren[segIndex].second];
    //Eq 19
    double bL = beta * sLeft.myBeta;
    double bR = beta * sRight.myBeta;
    updateRadius(sLeft.myIndex, bL);
    updateRadius(sRight.myIndex, bR);
    myVectSegments[myVectChildren[segIndex].first].myRadius = myVectSegments[1].myRadius*bL;
    myVectSegments[myVectChildren[segIndex].second].myRadius = myVectSegments[1].myRadius*bR;
  }
}


void
CoronaryArteryTree::updateRootRadius()
{
  unsigned int idRoot = 1;
  //Eq 18
  double Qpk = myVectSegments[idRoot].myKTerm*my_qTerm;
  double rRoot = pow((myVectSegments[idRoot].myResistance*Qpk)/my_pDrop,0.25);
  myVectSegments[idRoot].myRadius = rRoot;
  //updateRadius(idRoot);
  updateRadius(idRoot, 1.0);
}
/*
void
CoronaryArteryTree::updateSegmentRadiusToRoot(unsigned int segIndex)
{
  unsigned int idRoot = 1;
  double rRoot= myVectSegments[idRoot].myRadius;
  std::vector<unsigned int> v = getPathToRoot(myVectSegments[segIndex]);
  double radius = 1;
  for(size_t id=v.size()-1; id>=2; id--) {
    radius *= rRoot*myVectSegments[v[id]].myBeta;
    myVectSegments[v[id]].myRadius = radius;
  }
}
*/
double
CoronaryArteryTree::getLengthSegment(unsigned int segIndex)
{
  if(segIndex==0) return 0; //there is no parent segment
  DGtal::Z2i::RealPoint p1 = myVectSegments[segIndex].myCoordinate;
  DGtal::Z2i::RealPoint p2 = myVectSegments[myVectParent[segIndex]].myCoordinate;
  return (p1-p2).norm()*myLengthFactor;
}
/*
double
CoronaryArteryTree::computeTreeVolume(double mu, double lambda)
{
  double volume = 0;
  for(size_t idSeg = 1; idSeg<myVectSegments.size(); idSeg++) {
    volume += pow(getLengthSegment(idSeg),mu)*pow(myVectSegments[idSeg].myRadius,lambda);
  }
  return volume;
}
*/
/*
void
CoronaryArteryTree::updateScale(double scale)
{
  for(size_t it=0; it<this->myVectSegments.size(); it++) {
    this->myVectSegments[it].myCoordinate[0] *= scale;
    this->myVectSegments[it].myCoordinate[1] *= scale;
    if(it!=0) //except root radius
      this->myVectSegments[it].myRadius *= scale;
    //this->myVectSegments[it].myLength *= scale;
  }
  this->myRsupp *= scale;
  this->my_rPerf = this->myRsupp;
  //this->myDThresold=sqrt((M_PI*this->myRsupp*this->myRsupp)/this->myKTerm);
}
*/
double
CoronaryArteryTree::computeTotalVolume(unsigned int segIndex)
{
  //Volumn of current segment : Eq 20
  //double v = myVectSegments[segIndex].myRadius*getLengthSegment(segIndex);
  //double r1 = myVectSegments[segIndex].myRadius;
  //double r2 = myVectSegments[1].myRadius*myVectSegments[segIndex].myBeta;
  double l = getLengthSegment(segIndex);
  double r = myVectSegments[segIndex].myRadius;
  double v = myVectSegments[segIndex].myRadius*myVectSegments[segIndex].myRadius*getLengthSegment(segIndex);
  //std::cout<<"segIndex:"<<segIndex << " has Vol="<<v<<std::endl;
  double vL, vR;
  if(myVectChildren[segIndex].first!=0 && myVectChildren[segIndex].second!=0) {
    vL = computeTotalVolume(myVectChildren[segIndex].first);
    vR = computeTotalVolume(myVectChildren[segIndex].second);
    v += vL + vR;
    //std::cout<<"vL ("<<myVectChildren[segIndex].first<<") "<<vL<<" + vR("<<myVectChildren[segIndex].second<<") and "<<vR<<"="<<v<<std::endl;
  }
  //std::cout<<"Vol final ="<<v<<std::endl;
  return v;
}
/*
bool
CoronaryArteryTree::isAddable(const Point2D &p, unsigned int segIndex, unsigned int nbIter, unsigned int nbNeibour)
{
  DGtal::Z2i::RealPoint pOpt;
  
  Point2D newCenter = findBarycenter(p, segIndex);
  // Creation of the left child
  Segment<Point2D> sNewLeft;
  sNewLeft.myCoordinate = myVectSegments[segIndex].myCoordinate;
  sNewLeft.myIndex = myVectSegments.size();
  sNewLeft.myRadius = myVectSegments[segIndex].myRadius;
  sNewLeft.myFlow = myVectSegments[segIndex].myFlow;
  sNewLeft.myKTerm = myVectSegments[segIndex].myKTerm;
  sNewLeft.myResistance = myVectSegments[segIndex].myResistance;
  // Creation of the right child
  Segment<Point2D> sNewRight;
  sNewRight.myCoordinate = p;
  sNewRight.myIndex = myVectSegments.size()+1;
  sNewRight.myRadius = myVectSegments[segIndex].myRadius;
  sNewRight.myFlow = my_qTerm;
  sNewRight.myKTerm = 1;
  // Cretation a copy of center segment
  Segment<Point2D> sCurrent = myVectSegments[segIndex];
  sCurrent.myCoordinate = newCenter;
  sCurrent.myFlow = myVectSegments[segIndex].myFlow + my_qTerm;
  sCurrent.myKTerm = sNewLeft.myKTerm + 1;
  // Cretation a copy of parent segment
  Segment<Point2D> sParent = myVectSegments[myVectParent[segIndex]];
  double r0 = sCurrent.myRadius, r1 = sCurrent.myRadius, r2 = sCurrent.myRadius;
  bool res1 = kamyiaOptimization(sParent.myCoordinate, sCurrent, sNewLeft, sNewRight, nbIter, pOpt, r0, r1, r2);
  bool res2;
  if(res1) {
    res2 = isIntersecting(p, pOpt, segIndex, nbNeibour);
    if(!res2) {
      //Update optimal values
      sNewLeft.myRadius = r1;
      sNewRight.myRadius = r2;
      
      //Add left segment
      myVectSegments.push_back(sNewLeft);
      myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0));
      myVectParent.push_back(segIndex);
      
      // Update for the new parent of the new left segment
      unsigned int leftGrandChildIndex = myVectChildren[myVectSegments[segIndex].myIndex].first;
      unsigned int rightGrandChildIndex = myVectChildren[myVectSegments[segIndex].myIndex].second;
      myVectParent[leftGrandChildIndex] = sNewLeft.myIndex;
      myVectParent[rightGrandChildIndex] = sNewLeft.myIndex;
      
      // Update of the child of the new left segment (sNewLeft)
      myVectChildren[sNewLeft.myIndex].first = leftGrandChildIndex;
      myVectChildren[sNewLeft.myIndex].second = rightGrandChildIndex;
      
      //Add right segment
      myVectSegments.push_back(sNewRight);
      myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0));
      myVectParent.push_back(segIndex);
      myVectTerminals.push_back(sNewRight.myIndex);
      
      // Update center segment
      myVectSegments[segIndex].myCoordinate = pOpt;
      myVectSegments[segIndex].myRadius = r0;
      myVectSegments[segIndex].myKTerm = sCurrent.myKTerm;
      myVectSegments[segIndex].myFlow = sCurrent.myFlow;
      //update childrens of center segment
      myVectChildren[segIndex].first = sNewLeft.myIndex;
      myVectChildren[segIndex].second = sNewRight.myIndex;
      
      // update parameters
      myKTerm++;
      
      // Update physilogique paramaters
      updateFlow(myVectParent[segIndex]);
      updateResistanceTerminal(sNewRight.myIndex);
      updateResistance(sNewLeft.myIndex);
      updateRootRadius();
    }
  }
  return res1 && !res2;
}
*/
bool
CoronaryArteryTree::isAddable(const Point2D &p, unsigned int segIndex, unsigned int nbIter, double tolerance, unsigned int nbNeibour)
{
  DGtal::Z2i::RealPoint pOpt;
  
  Point2D newCenter = findBarycenter(p, segIndex);
  // Creation of the left child
  Segment<Point2D> sNewLeft;
  sNewLeft.myCoordinate = myVectSegments[segIndex].myCoordinate;
  sNewLeft.myIndex = myVectSegments.size();
  sNewLeft.myRadius = myVectSegments[segIndex].myRadius;
  sNewLeft.myFlow = myVectSegments[segIndex].myFlow;
  sNewLeft.myKTerm = myVectSegments[segIndex].myKTerm;
  sNewLeft.myResistance = myVectSegments[segIndex].myResistance;
  // Creation of the right child
  Segment<Point2D> sNewRight;
  sNewRight.myCoordinate = p;
  sNewRight.myIndex = myVectSegments.size()+1;
  sNewRight.myRadius = myVectSegments[segIndex].myRadius;
  sNewRight.myFlow = my_qTerm;
  sNewRight.myKTerm = 1;
  // Cretation a copy of center segment
  Segment<Point2D> sCurrent = myVectSegments[segIndex];
  sCurrent.myCoordinate = newCenter;
  sCurrent.myFlow = myVectSegments[segIndex].myFlow + my_qTerm;
  sCurrent.myKTerm = sNewLeft.myKTerm + 1;
  // Cretation a copy of parent segment
  Segment<Point2D> sParent = myVectSegments[myVectParent[segIndex]];
  DGtal::Z2i::RealPoint pCurrent = (sParent.myCoordinate+sNewLeft.myCoordinate)/2.0;
  double r0 = sCurrent.myRadius, r1 = sCurrent.myRadius, r2 = sCurrent.myRadius;
  bool res1 = true, res2 = false, isDone = false;
  double vol, volCurr, diffVol;
  CoronaryArteryTree cTreeCurr = *this;
  //std::cout<<"---------- segIndex: "<<segIndex<<std::endl;
  size_t i=0;
  for(size_t i=0; i<nbIter && res1 && !res2 && !isDone; i++) {
    res1 = kamyiaOptimization(pCurrent, sParent.myCoordinate, sCurrent, sNewLeft, sNewRight, 1, pOpt, r0, r1, r2);
    if(!res1) {
      if(volCurr>0) {
        isDone = true;
        *this = cTreeCurr;
      }
    }
    else {
      res2 = isIntersecting(p, pOpt, segIndex, nbNeibour, 2*cTreeCurr.myVectSegments[segIndex].myRadius);
      if(!res2) {
        CoronaryArteryTree cTree1 = *this;
        //Update optimal values
        sNewLeft.myRadius = r1;
        sNewRight.myRadius = r2;
        
        //Add left segment
        cTree1.myVectSegments.push_back(sNewLeft);
        cTree1.myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0));
        cTree1.myVectParent.push_back(segIndex);
        
        // Update for the new parent of the new left segment
        unsigned int leftGrandChildIndex = cTree1.myVectChildren[myVectSegments[segIndex].myIndex].first;
        unsigned int rightGrandChildIndex = cTree1.myVectChildren[myVectSegments[segIndex].myIndex].second;
        cTree1.myVectParent[leftGrandChildIndex] = sNewLeft.myIndex;
        cTree1.myVectParent[rightGrandChildIndex] = sNewLeft.myIndex;
        
        // Update of the child of the new left segment (sNewLeft)
        cTree1.myVectChildren[sNewLeft.myIndex].first = leftGrandChildIndex;
        cTree1.myVectChildren[sNewLeft.myIndex].second = rightGrandChildIndex;
        
        //Add right segment
        cTree1.myVectSegments.push_back(sNewRight);
        cTree1.myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0));
        cTree1.myVectParent.push_back(segIndex);
        cTree1.myVectTerminals.push_back(sNewRight.myIndex);
        
        // Update center segment
        cTree1.myVectSegments[segIndex].myCoordinate = pOpt;
        cTree1.myVectSegments[segIndex].myRadius = r0;
        cTree1.myVectSegments[segIndex].myKTerm = sCurrent.myKTerm;
        cTree1.myVectSegments[segIndex].myFlow = sCurrent.myFlow;
        //update childrens of center segment
        cTree1.myVectChildren[segIndex].first = sNewLeft.myIndex;
        cTree1.myVectChildren[segIndex].second = sNewRight.myIndex;
        
        // update parameters
        cTree1.myKTerm++;
        // Update physilogique paramaters
        cTree1.updateFlow(cTree1.myVectParent[segIndex]);
        cTree1.updateResistanceFromRoot();
        cTree1.updateRootRadius();
        
        vol = cTree1.computeTotalVolume();
        //std::cout<<"Iter "<<i<<" has tree Volume: "<< vol <<std::endl;
        if(i==0) {
          volCurr = vol;
          cTreeCurr = cTree1;
        }
        else {
          diffVol = volCurr - vol;
          if(fabs(diffVol) < tolerance) {
            //std::cout<<"Best volume at "<<i<<std::endl;
            isDone = true;
            *this = cTree1;
          }
          else {
            volCurr = vol;
            cTreeCurr = cTree1;
          }
        }
        //update new position, radius and flow
        pCurrent = pOpt;
        sCurrent.myRadius = cTree1.myVectSegments[segIndex].myRadius;
        sNewLeft.myRadius = cTree1.myVectSegments[cTree1.myVectChildren[segIndex].first].myRadius;
        sNewLeft.myFlow = cTree1.myVectSegments[cTree1.myVectChildren[segIndex].first].myFlow;
        sNewRight.myRadius = cTree1.myVectSegments[cTree1.myVectChildren[segIndex].second].myRadius;
        sNewRight.myFlow = cTree1.myVectSegments[cTree1.myVectChildren[segIndex].second].myFlow;
      }
    }
  }
  
  return res1 && !res2 && isDone;
}

bool
CoronaryArteryTree::addSegmentFromPoint(const Point2D &p,
                                        unsigned int nearIndex,
                                        double rLeft, double rRight)
{
  
  // To add a new segment we need to add two new points (p and middle point) with two segments.
  // (a) First point added : point associated to the nearest segment. (basic solution the center of the nearest)
  // (b) Second point added: the new point with its new segment.
  // to process (a): s
  Point2D newCenter = findBarycenter(p, nearIndex);
  
  // Creation of the left child
  Segment<Point2D> sNewLeft;
  sNewLeft.myCoordinate = myVectSegments[nearIndex].myCoordinate;
  sNewLeft.myRadius = rLeft;
  sNewLeft.myIndex = myVectSegments.size();
  sNewLeft.myFlow = myVectSegments[nearIndex].myFlow;
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
  Segment<Point2D> sNewRight;
  sNewRight.myCoordinate = p;
  sNewRight.myRadius = rRight;
  sNewRight.myIndex = myVectSegments.size();
  sNewRight.myFlow = my_qTerm;
  myVectSegments.push_back(sNewRight);
  myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0));
  myVectParent.push_back(nearIndex);
  myVectTerminals.push_back(sNewRight.myIndex);
  
  // Update center segment
  myVectSegments[nearIndex].myCoordinate = newCenter;
  myVectSegments[nearIndex].myFlow = sNewLeft.myFlow + sNewRight.myFlow;
  myVectSegments[nearIndex].myKTerm += 1;
  //update childrens of center segment
  myVectChildren[nearIndex].first = sNewLeft.myIndex;
  myVectChildren[nearIndex].second = sNewRight.myIndex;
  
  // Optimization
  bool resKamiya = kamyiaOptimization(nearIndex);
  
  // Update Kterm
  myKTerm++;
  
  // Update physilogique paramaters
  updateResistanceTerminal(sNewRight.myIndex);
  updateResistance(sNewLeft.myIndex);
  // Update root radius
  updateRootRadius();
  updateLengthFactor();
  
  return resKamiya;
}

bool
CoronaryArteryTree::isIntersecting(const Point2D &pNew, const Point2D &pCenter, unsigned int nearIndex, unsigned int nbNeibour, double minDistance)
{
  //bool inter = hasNearestIntersections(p, newCenter, 10);
  //bool inter = hasNearestIntersections(myVectParent[nearIndex], nearIndex, p, newCenter,  nbNeibour);
  bool inter = hasNearestIntersections(myVectParent[nearIndex], nearIndex, pNew, pCenter,  nbNeibour);
  if (inter){
    //DGtal::trace.warning() << "detection intersection" << std::endl;
    return true;
  }
  // Check barycenter is not too close considered segment
  if ((pCenter-myVectSegments[nearIndex].myCoordinate).norm() < minDistance ||
      (pCenter-myVectSegments[myVectParent[nearIndex]].myCoordinate).norm() < minDistance ||
      (pCenter - pNew).norm() < minDistance){
    //DGtal::trace.warning() << "new barycenter too close!!!!!!!!" << std::endl;
    return true;
  }
  // Check new point with new segment of barycenter:
  // - newPt and new segment [Barycenter-OriginNewSeg]
  // - newPt and new segment [Barycenter-FatherNewSeg]
  if ((pCenter-myVectSegments[nearIndex].myCoordinate).norm() < minDistance ||
      (pCenter-myVectSegments[myVectParent[nearIndex]].myCoordinate).norm() < minDistance) {
    //DGtal::trace.warning() << "initial too close to new!!!!!!!!" << std::endl;
    return true;
  }
  if (getProjDistance(nearIndex, pNew) < minDistance) {
    //DGtal::trace.warning() << "initial too close!!!!!!!!" << std::endl;
    return true;

  }
  if (isToCloseFromNearest(pNew, minDistance)){
    //DGtal::trace.warning() << "detection near too close " << std::endl;
    return true;
  }
  if (getProjDistance(nearIndex, pNew) < minDistance||
      getProjDistance(myVectParent[nearIndex], pNew) < minDistance||
      getProjDistance(myVectChildren[nearIndex].first, pNew) < minDistance||
      getProjDistance(myVectChildren[nearIndex].second, pNew) < minDistance){
    //DGtal::trace.warning() << "detection too close existing" << std::endl;
    return true;
  }
  return false;
}

// Test Version to debug...
/*
bool
CoronaryArteryTree::addSegmentFromPointBK(const Point2D &p,
                                        unsigned int nearIndex,
                                        double rLeft, double rRight)
{

  // To add a new segment we need to add two new points (p and middle point) with two segments.
  // (a) First point added : point associated to the nearest segment. (basic solution the center of the nearest)
  // (b) Second point added: the new point with its new segment.
  // to process (a): s
  Point2D newCenter = findBarycenter(p, nearIndex);
  //bool inter = hasNearestIntersections(p, newCenter, 10);
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
  Segment<Point2D> sNewLeft;
  sNewLeft.myCoordinate = myVectSegments[nearIndex].myCoordinate;
  sNewLeft.myRadius = rLeft;
  sNewLeft.myIndex = myVectSegments.size();
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
  Segment<Point2D> sNewRight;
  sNewRight.myCoordinate = p;
  sNewRight.myRadius = rRight;
  sNewRight.myIndex = myVectSegments.size();
  myVectSegments.push_back(sNewRight);
  myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0));
  myVectParent.push_back(nearIndex);
  myVectTerminals.push_back(sNewRight.myIndex);
  // Update center segment
  myVectSegments[nearIndex].myCoordinate = newCenter;
  //update childrens of center segment
  myVectChildren[nearIndex].first = sNewLeft.myIndex;
  myVectChildren[nearIndex].second = sNewRight.myIndex;

  myKTerm++;

  return true;
}
*/


void
CoronaryArteryTree::boardDisplay(double thickness, bool clearDisplay)
{
  if (clearDisplay){
    myBoard.clear();
  }
  // 57.5 from myBoard change scale
  double scaleBoard = 57.5;
  // drawing base circle
  //std::cout <<"My rsupp is " << myRsupp<<std::endl;
  myBoard.setPenColor(DGtal::Color::Blue);
  myBoard.setLineWidth(myVectSegments[0].myRadius*scaleBoard*thickness);
  myBoard.drawCircle(myTreeCenter[0], myTreeCenter[1], my_rPerf, 1);
  
  
  // draw root: (first segment is special reduced to one point and no parent).
  Point2D p0 = myVectSegments[0].myCoordinate;
  myBoard.setPenColor(DGtal::Color(10, 100, 0, 180));
  myBoard.setLineWidth(1.0);
  myBoard.fillCircle(p0[0], p0[1], myVectSegments[0].myRadius*thickness, 1);
  
  Point2D p1 = myVectSegments[1].myCoordinate;
  
  myBoard.setPenColor(DGtal::Color(180, 0, 0, 180));
  myBoard.setLineWidth(1.0);
  myBoard.fillCircle(p1[0], p1[1], myVectSegments[1].myRadius*thickness, 1);
  
  myBoard.setLineWidth(myVectSegments[1].myRadius*scaleBoard*thickness);
  myBoard.setPenColor(DGtal::Color(150, 0, 0, 150));
  myBoard.drawLine(p0[0], p0[1], p1[0], p1[1], 2);
  
  DGtal::GradientColorMap<int> cmap_grad( 0, myVectSegments.size()+1 );
  cmap_grad.addColor( DGtal::Color( 50, 50, 255 ) );
  cmap_grad.addColor( DGtal::Color( 255, 0, 0 ) );
  cmap_grad.addColor( DGtal::Color( 255, 255, 10 ) );
  unsigned int i = 0;
  for (auto s : myVectSegments)
  {
    // test if the segment is the root or its parent we do not display (already done).
    if (s.myIndex == 0 || s.myIndex == 1)
    {
      continue;
    }
    // distal node:
    Point2D distal = s.myCoordinate;
    Point2D proxital = myVectSegments[myVectParent[s.myIndex]].myCoordinate;
    myBoard.setLineWidth(myVectSegments[s.myIndex].myRadius*scaleBoard*thickness);
    //   myBoard.setPenColor(cmap_grad(i));
    myBoard.setPenColor(DGtal::Color::Gray);
    myBoard.drawLine(distal[0], distal[1], proxital[0], proxital[1],2);
    //      std::cout << " myRsupp Value = "<< myRsupp<<std::endl;
    myBoard.setLineWidth(1.0);
    myBoard.setPenColor(DGtal::Color(180, 0, 0, 180));
    myBoard.fillCircle(distal[0], distal[1], myVectSegments[s.myIndex].myRadius*thickness, 1 );
    i++;
    
  }
}


void
CoronaryArteryTree::exportBoardDisplay(const std::string &fileName,
                                       double thickness,
                                       bool updateDisplay,
                                       bool clearDisplay){
  if (updateDisplay){
    boardDisplay(thickness, clearDisplay);
  }
  std::string ext = fileName.substr(fileName.find_last_of(".") + 1);
  if (ext == "svg")
  {
    myBoard.saveSVG(fileName.c_str());
  }
  else if (ext == "eps") {
    myBoard.saveEPS(fileName.c_str());
  }
}

/*
CoronaryArteryTree::Point2D
CoronaryArteryTree::getSegmentCenter(unsigned int i)
{
  return getSegmentCenter(myVectSegments.at(i));
}
*/
/*
CoronaryArteryTree::Point2D
CoronaryArteryTree::getSegmentCenter(const Segment<Point2D> &s)
{
  Point2D res;
  res=s.myCoordinate;
  // Test if the segment is not a the special root segment (reducted to one point)
  if (s.myIndex != 0)
  {
    unsigned int p = myVectParent.at(s.myIndex);
    res += myVectSegments.at(p).myCoordinate;
    res[0] /= 2.0;
    res[1] /= 2.0;
  }
  return res;
}
*/

unsigned int
CoronaryArteryTree::getNearestSegment(const Point2D &pt)
{
  unsigned int sNear=1;
  double distMin = getProjDistance(myVectSegments[0].myCoordinate,
                                   myVectSegments[1].myCoordinate, pt);
  for (const auto &s: myVectSegments)
  {
    if (s.myIndex==0)
      continue;
    
    double d = getProjDistance(s.myCoordinate,
                              myVectSegments[myVectParent[s.myIndex]].myCoordinate, pt);

    if (d < distMin)
    {
      distMin = d;
      sNear = s.myIndex;
    }
  }
  return sNear;
}


/*
unsigned int
CoronaryArteryTree::getParentSegment(const Segment<Point2D> &s){
  return myVectParent[s.myIndex];
}
*/
/*
unsigned int
CoronaryArteryTree::getLeftChild(const Segment<Point2D> &s){
  return myVectChildren[s.myIndex].first;
}
*/
/*
unsigned int
CoronaryArteryTree::getRightChild(const Segment<Point2D> &s){
  return myVectChildren[s.myIndex].second;
}
*/
/*
double
CoronaryArteryTree::compDistCriteria(const Point2D &p, unsigned int indexNode)
{
  Point2D pDist = myVectSegments[indexNode].myCoordinate;
  Point2D pProxi = myVectSegments[myVectParent[indexNode]].myCoordinate;
  Point2D pProj;
  bool isInside = projectOnStraightLine(pDist, pProxi, p, pProj);
  if (isInside)
  {
    return (pProj-p).norm();
  }
  return std::min((p-pDist).norm(), (p-pProxi).norm());
}
*/
std::vector<unsigned int>
CoronaryArteryTree::getPathToRoot(const Segment<Point2D> &s)
{
  std::vector<unsigned int> res;
  Segment<Point2D> sC = s;
  while (sC.myIndex != 0){
    res.push_back(sC.myIndex);
    sC = myVectSegments[ myVectParent[sC.myIndex]];
  }
  return res;
  
}

/*
void
CoronaryArteryTree::udpatePerfusionArea(){
  myCurrAPerf = (myKTerm+1.0)*M_PI*myRsupp*myRsupp;
  
}
*/

CoronaryArteryTree::Point2D
CoronaryArteryTree::generateNewLocation(unsigned int nbTrials){
  Point2D res;
  double myDThresold = getDistanceThreshold();
  bool found = false;
  unsigned int n = nbTrials;
  while(!found && n >= 0) {
    n--;
    auto p = generateALocation(myDThresold);
    found = p.second;
    if (found) {
      res = p.first;
    }
    if (n==0) {
      n = nbTrials;
      myDThresold *= 0.9;
    }
  }
  
  return res;
}


std::pair<CoronaryArteryTree::Point2D, bool>
CoronaryArteryTree::generateALocation(double myDThresold) {
  Point2D res = generateRandomPtOnDisk(myTreeCenter, my_rPerf);
  bool isComp = true;
  unsigned int id = 1;
  while ( isComp && id < myVectTerminals.size() ) {
    isComp = (myVectSegments[myVectTerminals[id]].myCoordinate - res).norm() > myDThresold;
    id++;
  }
  return  std::pair<Point2D, bool> {res, isComp};
}
/*
double
CoronaryArteryTree::getDistance(unsigned int index,
                                const Point2D &p ) const {
  return (myVectSegments[index].myCoordinate-p).norm();
}
*/
bool
CoronaryArteryTree::isToCloseFromNearest(const Point2D &p, double minDist) const{
  double d = getProjDistance(getN_NearestSegments(p,1)[0],p);
  return d < minDist;
}
double
CoronaryArteryTree::getProjDistance(const Point2D &p0, const Point2D &p1, const Point2D &p) const{
 // return ((p0+p1)/2.0 - p).norm();
  double result = 0.0;
  Point2D pProj(0,0);
  bool isInside = projectOnStraightLine(p0, p1, p, pProj);
  if (isInside)
  {
    result = (pProj-p).norm();
  }
  else
  {
    result = std::min((p0-p).norm(), (p1-p).norm());
  }

  return result;
}

double
CoronaryArteryTree::getProjDistance(unsigned int index, const Point2D &p ) const {
  Point2D p0 = myVectSegments[index].myCoordinate;
  Point2D p1 = myVectSegments[myVectParent[index]].myCoordinate;
  return  getProjDistance(p0, p1, p);
}



std::vector<unsigned int>
CoronaryArteryTree::getN_NearestSegments(const Point2D &p, unsigned int n) const {
  std::vector<unsigned int> res;
  res.push_back(1);
  for (unsigned int i=2; i < myVectSegments.size(); i++){
    double d = getProjDistance(i, p);
    std::vector<unsigned int>::iterator it = res.end();
    int k = (int)(res.size())-1;
    while ( k >= 0 && getProjDistance(res[k], p)>=d ){
      k--;
      it--;
    }
    res.insert(it, i);
    if (res.size() > n){
      res.pop_back();
    }
  }
  return res;
}


bool
CoronaryArteryTree::hasNearestIntersections(const Point2D &p0,
                                            const Point2D &p1, unsigned int n) const {
  Point2D b = (p1+p0)/2;
  std::vector<unsigned int> near = getN_NearestSegments(b, n);
  for (const auto &s : near){
    if (hasIntersection(p0, p1, myVectSegments[s].myCoordinate,
                        myVectSegments[myVectParent[s]].myCoordinate ))
    {
      return  true;
    }
  }
  return false;
}
bool
CoronaryArteryTree::hasNearestIntersections(unsigned int indexPFather,
                                            unsigned int indexPChild,
                                            const Point2D &pAdded,
                                            const Point2D &pBifurcation, unsigned int n) const{
  std::vector<unsigned int> near = getN_NearestSegments(pAdded, n);
  Point2D p0  = myVectSegments[indexPFather].myCoordinate;
  Point2D p1  = myVectSegments[indexPChild].myCoordinate;
  
  for (const auto &s : near){
    // ignoring self initial segment.
    if (myVectParent[s] == indexPFather && s == indexPChild){
      continue;
    }
    bool inter1 = hasIntersection(p0, pBifurcation, myVectSegments[s].myCoordinate,
                                  myVectSegments[myVectParent[s]].myCoordinate  );
    bool inter2 = hasIntersection(p1, pBifurcation, myVectSegments[s].myCoordinate,
                                  myVectSegments[myVectParent[s]].myCoordinate  );
    bool inter3 = hasIntersection(pAdded, pBifurcation, myVectSegments[s].myCoordinate,
                                  myVectSegments[myVectParent[s]].myCoordinate  );

    if (inter1 || inter2 || inter3)
    {
      return  true;
    }
  }
  return false;
  
}



bool
CoronaryArteryTree::kamyiaOptimization(unsigned int index,
                                       unsigned int nbIter){
  
  // Prepare
  Segment<Point2D> sParent = myVectSegments[myVectParent[index]];
  Segment<Point2D> sL = myVectSegments[myVectChildren[index].first];
  Segment<Point2D> sR = myVectSegments[myVectChildren[index].second];
  
  DGtal::Z2i::RealPoint pParent = sParent.myCoordinate;
  DGtal::Z2i::RealPoint pL = sL.myCoordinate;
  DGtal::Z2i::RealPoint pR = sR.myCoordinate;
  
  //std::cout<<"myFlow="<< myVectSegments[index].myFlow<<" and sL.myFlow="<<sL.myFlow<<" and sR.myFlow="<<sR.myFlow<<std::endl;
  double ratioQ = sL.myFlow/sParent.myFlow;//0.5;
  std::cout<<"ratioQ="<< sL.myFlow/myVectSegments[index].myFlow<<std::endl;
  double r0 = myVectSegments[index].myRadius;
  double R0 = r0*r0; //R0 = r0*r0
  double R1 = r0*r0; //R1 = r1*r1
  double R2 = r0*r0; //R2 = r2*r2
  
  double f1 = sL.myFlow;//ratioQ*f0;//k * r1*r1*r1; //left
  double f2 = sR.myFlow;//(1.0-ratioQ)*f0;//k * r2*r2*r2; //right
  double f0 =  f1 + f2;//R0*R0*R0; //middle
  // Starting position from Equation (21) for initialisation as mentionned page 11 [Clara Jaquet et HT]
  //DGtal::Z2i::RealPoint pb ((f0*pParent[0]+f1*pL[0]+f2*pR[0])/(2.0*f0), (f0*pParent[1]+f1*pL[1]+f2*pR[1])/(2.0*f0));
  DGtal::Z2i::RealPoint pb ((pParent[0]+pL[0])/(2.0), (pParent[1]+pL[1])/(2.0)); //center of the old segment
  //std::cout<<"Init:"<<pb<<std::endl;
  double l0 = (pParent - pb).norm()*myLengthFactor;
  double l1 = (pL - pb).norm()*myLengthFactor;
  double l2 = (pR - pb).norm()*myLengthFactor;
  //std::cout<<"l0="<<l0<<", l1="<<l1<<", l2="<<l2<<std::endl;
  //double kappa = 8 * my_nu / M_PI;
  //double deltaP1k = kappa*((f0*l0)/(R0*R0)+(f1*l1)/(R1*R1));
  //double deltaP2k = kappa*((f0*l0)/(R0*R0)+(f2*l2)/(R2*R2));
  double deltaP1 = (f0*l0)/(R0*R0)+(f1*l1)/(R1*R1);
  double deltaP2 = (f0*l0)/(R0*R0)+(f2*l2)/(R2*R2);
  
  double rr1 = R1;
  double rr2 = R2;
  // The variable to solve for with its initial value. It will be
  // mutated in place by the solver.
  //DGtal::Z2i::RealPoint pbInit ((f0*pParent[0]+f1*pL[0]+f2*pR[0])/(2.0*f0), (f0*pParent[1]+f1*pL[1]+f2*pR[1])/(2.0*f0));
  //DGtal::trace.info() << pbInit << std::endl;
  
  bool hasSolution = true;
  
  for (int i = 0; i<nbIter && hasSolution; i++)
  {
    //DGtal::trace.progressBar(i, nbIter);
    hasSolution = kamiyaOpt(my_gamma, deltaP1, deltaP2, f0, f1, f2, l0, l1, l2, rr1, rr2);
    // Equation 27
    R0 = pow(f0*(pow(rr1, my_gamma)/f1 + pow(rr2, my_gamma)/f2), 1.0/my_gamma);
    R1 = rr1;
    R2 = rr2;
    // Equation (26) page 13
    pb[0] = (pParent[0]*R0/l0 + pL[0]*R1/l1 + pR[0]*R2/l2)/(R0/l0+R1/l1+R2/l2);
    pb[1] = (pParent[1]*R0/l0 + pL[1]*R1/l1 + pR[1]*R2/l2)/(R0/l0+R1/l1+R2/l2);
    //Update values
    deltaP1 = (f0*l0)/(R0*R0)+(f1*l1)/(R1*R1);
    deltaP2 = (f0*l0)/(R0*R0)+(f2*l2)/(R2*R2);
    l0 = (pParent - pb).norm()*myLengthFactor;
    l1 = (pL - pb).norm()*myLengthFactor;
    l2 = (pR - pb).norm()*myLengthFactor;
    std::cout << "r0 : " << sqrt(R0) << ", r1 : " << sqrt(R1) << ", sqrt(R2) : " << sqrt(R2) << "\n";
    std::cout << "pbNew[0]  : " << pb[0] << " and pbNew[1] " << pb[1] << "\n";
  }
  // If there is a solution, then update result
  if(hasSolution) {
    // Update center point coordinate
    myVectSegments[index].myCoordinate[0] = pb[0];
    myVectSegments[index].myCoordinate[1] = pb[1];
    double mR0 = sqrt(R0);
    double mR1 = sqrt(R1);
    double mR2 = sqrt(R2);
    double mL0 = l0;
    double mL1 = l1;
    double mL2 = l2;
    
    // Update radius of center, left and right segments
    myVectSegments[index].myRadius = mR0;
    myVectSegments[myVectChildren[index].first].myRadius =  mR1;
    myVectSegments[myVectChildren[index].second].myRadius =  mR2;
    
    if(2*mR0<=mL0 && 2*mR1<=mL1 && 2*mR2<mL2)
      hasSolution = true;
    else
      hasSolution = false;
  }
  return hasSolution;
}

bool
CoronaryArteryTree::kamyiaOptimization(const DGtal::Z2i::RealPoint& pParent,
                                       const Segment<Point2D>& sCurrent,
                                       const Segment<Point2D>& sL,
                                       const Segment<Point2D>& sR,
                                       unsigned int nbIter,
                                       DGtal::Z2i::RealPoint& pOpt,
                                       double& r0, double& r1, double& r2){
  
  // Prepare
  DGtal::Z2i::RealPoint pL = sL.myCoordinate;
  DGtal::Z2i::RealPoint pR = sR.myCoordinate;
  
  //std::cout<<"myFlow="<< myVectSegments[index].myFlow<<" and sL.myFlow="<<sL.myFlow<<" and sR.myFlow="<<sR.myFlow<<std::endl;
  double ratioQ = sL.myFlow/sCurrent.myFlow;//0.5;
  //std::cout<<"ratioQ="<< ratioQ<<std::endl;
  double rr0 = sCurrent.myRadius;
  double R0 = rr0*rr0; //R0 = r0*r0
  double R1 = rr0*rr0; //R1 = r1*r1
  double R2 = rr0*rr0; //R2 = r2*r2
  
  double f1 = sL.myFlow;//ratioQ*f0;//k * r1*r1*r1; //left
  double f2 = sR.myFlow;//(1.0-ratioQ)*f0;//k * r2*r2*r2; //right
  double f0 =  f1 + f2;//R0*R0*R0; //middle
  // Starting position from Equation (21) for initialisation as mentionned page 11 [Clara Jaquet et HT]
  //DGtal::Z2i::RealPoint pb ((f0*pParent[0]+f1*pL[0]+f2*pR[0])/(2.0*f0), (f0*pParent[1]+f1*pL[1]+f2*pR[1])/(2.0*f0));
  DGtal::Z2i::RealPoint pb ((pParent[0]+pL[0])/(2.0), (pParent[1]+pL[1])/(2.0)); //center of the old segment
  //std::cout<<"Init:"<<pb<<std::endl;
  double l0 = (pParent - pb).norm()*myLengthFactor;
  double l1 = (pL - pb).norm()*myLengthFactor;
  double l2 = (pR - pb).norm()*myLengthFactor;
  //std::cout<<"l0="<<l0<<", l1="<<l1<<", l2="<<l2<<std::endl;
  //double kappa = 8 * my_nu / M_PI;
  //double deltaP1k = kappa*((f0*l0)/(R0*R0)+(f1*l1)/(R1*R1));
  //double deltaP2k = kappa*((f0*l0)/(R0*R0)+(f2*l2)/(R2*R2));
  double deltaP1 = (f0*l0)/(R0*R0)+(f1*l1)/(R1*R1);
  double deltaP2 = (f0*l0)/(R0*R0)+(f2*l2)/(R2*R2);
  
  double rr1 = R1;
  double rr2 = R2;
  
  // The variable to solve for with its initial value. It will be
  // mutated in place by the solver.
  
  bool hasSolution = true;
  
  for (int i = 0; i<nbIter && hasSolution; i++) {
    //DGtal::trace.progressBar(i, nbIter);
    hasSolution = kamiyaOpt(my_gamma, deltaP1, deltaP2, f0, f1, f2, l0, l1, l2, rr1, rr2);
    // Equation 27
    R0 = pow(f0*(pow(rr1, my_gamma)/f1 + pow(rr2, my_gamma)/f2), 1.0/my_gamma);
    R1 = rr1;
    R2 = rr2;
    // Equation (26) page 13
    pb[0] = (pParent[0]*R0/l0 + pL[0]*R1/l1 + pR[0]*R2/l2)/(R0/l0+R1/l1+R2/l2);
    pb[1] = (pParent[1]*R0/l0 + pL[1]*R1/l1 + pR[1]*R2/l2)/(R0/l0+R1/l1+R2/l2);
    //Update values
    deltaP1 = (f0*l0)/(R0*R0)+(f1*l1)/(R1*R1);
    deltaP2 = (f0*l0)/(R0*R0)+(f2*l2)/(R2*R2);
    l0 = (pParent - pb).norm()*myLengthFactor;
    l1 = (pL - pb).norm()*myLengthFactor;
    l2 = (pR - pb).norm()*myLengthFactor;
    //std::cout << "r0 : " << sqrt(R0) << ", r1 : " << sqrt(R1) << ", sqrt(R2) : " << sqrt(R2) << "\n";
    //std::cout << "R0 : " << R0 << ", R1 : " << R1 << ", R2 : " << R2 << "\n";
    //std::cout << "pbNew[0]  : " << pb[0] << " and pbNew[1] " << pb[1] << "\n";
  }
  
  // If there is a solution, then update result
  if(hasSolution) {
    pOpt = pb;
    r0 = sqrt(R0);
    r1 = sqrt(R1);
    r2 = sqrt(R2);
    double mL0 = l0;//(pb - pParent).norm();
    double mL1 = l1;//(pb - pL).norm();
    double mL2 = l2;//(pb - pR).norm();
    
    if(2*r0<=mL0 && 2*r1<=mL1 && 2*r2<mL2)
      hasSolution = true;
    else
      hasSolution = false;
  }
  return hasSolution;
}

bool
CoronaryArteryTree::kamyiaOptimization(const DGtal::Z2i::RealPoint& pCurrent,
                                       const DGtal::Z2i::RealPoint& pParent,
                                       const Segment<Point2D>& sCurrent,
                                       const Segment<Point2D>& sL,
                                       const Segment<Point2D>& sR,
                                       unsigned int nbIter,
                                       DGtal::Z2i::RealPoint& pOpt,
                                       double& r0, double& r1, double& r2){
  
  // Prepare
  DGtal::Z2i::RealPoint pL = sL.myCoordinate;
  DGtal::Z2i::RealPoint pR = sR.myCoordinate;
  
  //std::cout<<"myFlow="<< myVectSegments[index].myFlow<<" and sL.myFlow="<<sL.myFlow<<" and sR.myFlow="<<sR.myFlow<<std::endl;
  //double ratioQ = sL.myFlow/sCurrent.myFlow;//0.5;
  //std::cout<<"ratioQ="<< ratioQ<<std::endl;
  //double rr0 = sCurrent.myRadius;
  double R0 = sCurrent.myRadius*sCurrent.myRadius;//rr0*rr0; //R0 = r0*r0
  double R1 = sL.myRadius*sL.myRadius;//rr0*rr0; //R1 = r1*r1
  double R2 = sR.myRadius*sR.myRadius;//rr0*rr0; //R2 = r2*r2
  
  double f1 = sL.myFlow;//ratioQ*f0;//k * r1*r1*r1; //left
  double f2 = sR.myFlow;//(1.0-ratioQ)*f0;//k * r2*r2*r2; //right
  double f0 =  f1 + f2;//R0*R0*R0; //middle
  // Starting position from Equation (21) for initialisation as mentionned page 11 [Clara Jaquet et HT]
  //DGtal::Z2i::RealPoint pb ((f0*pParent[0]+f1*pL[0]+f2*pR[0])/(2.0*f0), (f0*pParent[1]+f1*pL[1]+f2*pR[1])/(2.0*f0));
  //DGtal::Z2i::RealPoint pb ((pParent[0]+pL[0])/(2.0), (pParent[1]+pL[1])/(2.0)); //center of the old segment
  DGtal::Z2i::RealPoint pb (pCurrent[0], pCurrent[1]);
  //std::cout<<"Init:"<<pb<<std::endl;
  double l0 = (pParent - pb).norm()*myLengthFactor;
  double l1 = (pL - pb).norm()*myLengthFactor;
  double l2 = (pR - pb).norm()*myLengthFactor;
  //std::cout<<"l0="<<l0<<", l1="<<l1<<", l2="<<l2<<std::endl;
  double kappa = 8 * my_nu / M_PI;
  double deltaP1k = kappa*((f0*l0)/(R0*R0)+(f1*l1)/(R1*R1));
  double deltaP2k = kappa*((f0*l0)/(R0*R0)+(f2*l2)/(R2*R2));
  double deltaP1 = (f0*l0)/(R0*R0)+(f1*l1)/(R1*R1);
  double deltaP2 = (f0*l0)/(R0*R0)+(f2*l2)/(R2*R2);
  
  double rr1 = R1;
  double rr2 = R2;
  
  // The variable to solve for with its initial value. It will be
  // mutated in place by the solver.
  
  bool hasSolution = true;
  
  for (int i = 0; i<nbIter && hasSolution; i++) {
    //DGtal::trace.progressBar(i, nbIter);
    hasSolution = kamiyaOpt(my_gamma, deltaP1, deltaP2, f0, f1, f2, l0, l1, l2,rr1, rr2);
    // Equation 27
    R0 = pow(f0*(pow(rr1, my_gamma)/f1 + pow(rr2, my_gamma)/f2), 1.0/my_gamma);
    R1 = rr1;
    R2 = rr2;
    // Equation (26) page 13
    pb[0] = (pParent[0]*R0/l0 + pL[0]*R1/l1 + pR[0]*R2/l2)/(R0/l0+R1/l1+R2/l2);
    pb[1] = (pParent[1]*R0/l0 + pL[1]*R1/l1 + pR[1]*R2/l2)/(R0/l0+R1/l1+R2/l2);
    //Update values
    deltaP1 = (f0*l0)/(R0*R0)+(f1*l1)/(R1*R1);
    deltaP2 = (f0*l0)/(R0*R0)+(f2*l2)/(R2*R2);
    l0 = (pParent - pb).norm()*myLengthFactor;
    l1 = (pL - pb).norm()*myLengthFactor;
    l2 = (pR - pb).norm()*myLengthFactor;
    //std::cout << "r0 : " << sqrt(R0) << ", r1 : " << sqrt(R1) << ", sqrt(R2) : " << sqrt(R2) << "\n";
    //std::cout << "R0 : " << R0 << ", R1 : " << R1 << ", R2 : " << R2 << "\n";
    //std::cout << "pbNew[0]  : " << pb[0] << " and pbNew[1] " << pb[1] << "\n";
  }
  
  // If there is a solution, then update result
  if(hasSolution) {
    pOpt = pb;
    r0 = sqrt(R0);
    r1 = sqrt(R1);
    r2 = sqrt(R2);
    //double mL0 = (pb - pParent).norm();
    //double mL1 = (pb - pL).norm();
    //double mL2 = (pb - pR).norm();
    
    //if(2*r0<=mL0 && 2*r1<=mL1 && 2*r2<mL2)
    if(2*r0<=l0 && 2*r1<=l1 && 2*r2<l2)
      hasSolution = true;
    else
      hasSolution = false;
  }
  return hasSolution;
}

void
CoronaryArteryTree::selfDisplay( std::ostream & out ) const {
  out << std::endl << "----" << std::endl;
  out << "CoronaryArteryTree: " << std::endl;
  out << "main parameters: myKTerm: "  << myKTerm << "\n\t myRsupp: " << myRsupp << std::endl;
  out << "\n\t Nb Segments from container: " << myVectSegments.size() << std::endl;
  out << "----" << std::endl;
}


std::ostream&
operator<< ( std::ostream & out,
            const CoronaryArteryTree & aCoronaryTree )
{
  aCoronaryTree.selfDisplay ( out );
  return out;
}

bool
CoronaryArteryTree::addSegmentFromPointWithBarycenter(const Point2D &p)
{
  unsigned int nearIndex = getNearestSegment(p);
  return addSegmentFromPointWithBarycenter(p, nearIndex);
}

bool
CoronaryArteryTree::addSegmentFromPointWithBarycenter(const Point2D &p, unsigned int nearIndex)
{
  // (a) First point added : point associated to the nearest segment. (basic solution the center of the nearest)
  Point2D newCenter = findBarycenter(p, nearIndex);
  
  // Creation of the left child
  Segment<Point2D> sNewLeft;
  sNewLeft.myCoordinate = myVectSegments[nearIndex].myCoordinate;
  sNewLeft.myRadius = 1.0;
  sNewLeft.myIndex = myVectSegments.size();
  sNewLeft.myFlow = myVectSegments[nearIndex].myFlow;
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
  Segment<Point2D> sNewRight;
  sNewRight.myCoordinate = p;
  sNewRight.myRadius = 1.0;
  sNewRight.myIndex = myVectSegments.size();
  sNewRight.myFlow = my_qTerm;
  myVectSegments.push_back(sNewRight);
  myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0));
  myVectParent.push_back(nearIndex);
  myVectTerminals.push_back(sNewRight.myIndex);
  
  // Update center segment
  myVectSegments[nearIndex].myCoordinate = newCenter;
  //update childrens of center segment
  myVectChildren[nearIndex].first = sNewLeft.myIndex;
  myVectChildren[nearIndex].second = sNewRight.myIndex;
  
  // Update Kterm
  myKTerm++;
  
  //updateRadius();
  updateRootRadius();
  return true;
  
}
/*
double
CoronaryArteryTree::dProjCalculation(const Point2D &p,unsigned int Index )
{
  // p is the point that we want to add.
  
  Point2D distalPoint=myVectSegments[Index].myCoordinate;
  Point2D proximalPoint = myVectSegments[myVectParent[Index]].myCoordinate;
  double length  = pow(pow((proximalPoint[0] - distalPoint[0]), 2.0) + pow((proximalPoint[1] - distalPoint[1]), 2.0),0.5);
  double dProj = (((distalPoint[0]-proximalPoint[0])*(p[0] - proximalPoint[0]))+((distalPoint[1]-proximalPoint[1])*(p[1] - proximalPoint[1])))/length;
  
  return dProj;
  
}
*/
/*
double
CoronaryArteryTree::dCritCalculation(const Point2D &p,unsigned int Index )
{
  Point2D distalPoint=myVectSegments[Index].myCoordinate;
  Point2D proximalPoint = myVectSegments[myVectParent[Index]].myCoordinate;
  double length  = pow(pow((proximalPoint[0] - distalPoint[0]), 2.0) + pow((proximalPoint[1] - distalPoint[1]), 2.0),0.5); // we had to add length to segment
  double dCrit;
  if(dProjCalculation(p,Index) >0 && dProjCalculation(p,Index)< 1)
  {
    dCrit = abs(((distalPoint[1] - proximalPoint[1])*(p[0]-proximalPoint[0]))+((proximalPoint[0]-distalPoint[0])*(p[1]-proximalPoint[1])))/length;
  }
  else
  {
    dCrit = std::min(pow(pow((p[0]-proximalPoint[0]),2)+ pow(p[1]-proximalPoint[1],2 ),0.5),  pow(pow((p[0]-distalPoint[0]),2)+ pow(p[1]-distalPoint[1],2 ),0.5) );
  }
  return dCrit;
  
}
*/
double
CoronaryArteryTree::getDistanceThreshold()
{
  //Eq 12
  return sqrt((M_PI*my_rPerf*my_rPerf)/myKTerm);
  //return sqrt((M_PI*myRsupp*myRsupp)/myKTerm);
}

/*
bool
CoronaryArteryTree::updateRadius()
{
  for (auto s : myVectSegments)
  {
    myVectSegments[s.myIndex].myRadius=sqrt(my_qTerm/(M_PI*GetLength(s.myIndex)));
  }
  //    cout << "max gene = " << MaxGene<<endl;
  return TRUE;
  
}
*/
/*
bool
CoronaryArteryTree::updateRadius2(unsigned int index )
{
  //First, we have to change the radius of the two new segments. They have to carry the flow qTerm;
  Segment<Point2D> parent;
  parent = myVectSegments[myVectParent[index]];
  
  
  myVectSegments[getLeftChild(parent)].myRadius = sqrt(1*my_qTerm/(M_PI*GetLength(getLeftChild(parent))));
  
  myVectSegments[getRightChild(parent)].myRadius = sqrt((1+1-1)*my_qTerm/(M_PI*GetLength(getRightChild(parent))));
  // Now, we need to change the parents according to the Bifurcation rule;
  
  
  myVectSegments[parent.myIndex].myRadius =sqrt((myKTerm)*my_qTerm/(M_PI*GetLength(parent.myIndex)));
  
  return TRUE;
  
}
*/





/*
double
CoronaryArteryTree::GetLength(unsigned int Index)
{
  Point2D distalPoint=myVectSegments[Index].myCoordinate;
  Point2D proximalPoint = myVectSegments[myVectParent[Index]].myCoordinate;
  double length  = pow(pow((proximalPoint[0] - distalPoint[0]), 2.0) + pow((proximalPoint[1] - distalPoint[1]), 2.0),0.5);
  return length;
}
*/

/*
double
CoronaryArteryTree::GetTotalVolume(const Point2D &p1,const Point2D &p2,const Point2D &p3,const Point2D &pOpti)
{
  double H1,H2,H3;
  
  H1 = sqrt(pow(p1[0]-pOpti[0],2)+pow(p1[1]-pOpti[1],2));
  H2 = sqrt(pow(p2[0]-pOpti[0],2)+pow(p2[1]-pOpti[1],2));
  H3 = sqrt(pow(p3[0]-pOpti[0],2)+pow(p3[1]-pOpti[1],2));
  
  return M_PI*H1+M_PI*H2+M_PI*H3;
  
}
*/
bool
CoronaryArteryTree::addSegment(const Point2D &NewPoint,const Point2D &OptimizePoint, unsigned int nearIndex)
{
  
  // To add a new segment we need to add two new points (p and middle point) with two segments.
  // (a) First point added : point associated to the nearest segment. (basic solution the center of the nearest)
  // (b) Second point added: the new point with its new segment.
  
  // to process (a): s middle
  Segment<Point2D> sMiddle;
  sMiddle.myCoordinate = OptimizePoint;
  sMiddle.myRadius = 0.1;
  sMiddle.myIndex = myVectSegments.size();
  myVectSegments.push_back(sMiddle);
  myVectParent.push_back(myVectParent[nearIndex]);
  myVectChildren.push_back(SegmentChildren(nearIndex,  myVectSegments.size()));
  
  // to process (b): new point to Optimal position
  Segment<Point2D> sNew;
  sNew.myCoordinate = NewPoint;
  sNew.myRadius = 0.1;
  sNew.myIndex = myVectSegments.size();
  myVectSegments.push_back(sNew);
  myVectParent.push_back(sMiddle.myIndex);
  // update parent with new (a).
  myVectParent[nearIndex] = sMiddle.myIndex;
  
  
  return true;
}

CoronaryArteryTree::Point2D
CoronaryArteryTree::findBarycenter(const Point2D &p, unsigned int index)
{
  Point2D first_point=p;
  Point2D second_point=myVectSegments[index].myCoordinate;
  Point2D third_point=myVectSegments[myVectParent[index]].myCoordinate;
  Point2D barycenter((first_point[0]+second_point[0]+third_point[0])/3.0,(first_point[1]+second_point[1]+third_point[1])/3.0);
  return barycenter;
}

/*
int
CoronaryArteryTree::AddFirstSegmentonImage()
{
  DGtal::Z2i::RealPoint p(320,200);
  
  Segment<Point2D> s;
  s.myCoordinate = p;
  s.myRadius = 0.1;
  s.myIndex = myVectSegments.size();
  myVectSegments.push_back(s);
  myVectParent.push_back(0);
  myVectChildren.push_back(SegmentChildren(0,0));
  return true;
}
*/
/*
double
CoronaryArteryTree::FindXmax(int xDim, int yDim)
{
  
  return myRsupp*cos(atan(xDim/yDim));
}
*/
/*
double
CoronaryArteryTree::FindYmax(int xDim, int yDim)
{
  
  return myRsupp*sin(atan(xDim/yDim));
  
}
*/
/*
// The following function converts a M_PIxelfrom the image into a point in the supporting circle
CoronaryArteryTree::Point2D
CoronaryArteryTree::fromImageToCircle(int ximage,int yimage,int xdim,int ydim)
{
  double x,y;
  
  double xmax,ymax;
  xmax = FindXmax(xdim,ydim);
  ymax = FindYmax(xdim,ydim);
  
  x = -ximage*xmax*2/xdim + xmax;
  y = -yimage*ymax*2/ydim + ymax;
  DGtal::Z2i::RealPoint NewPoint(x,y);
  
  return NewPoint;
  
}
*/

/*
// The following function converts a random generated point into a couple (x,y) in the M_PIcture
CoronaryArteryTree::Point2D
CoronaryArteryTree::fromCircleToImage(string fileName, double x, double y, int xdim,int ydim )
{
  
  
  double xmax,ymax;
  xmax = FindXmax(xdim,ydim);
  ymax = FindYmax(xdim,ydim);
  
  
  int ximage = -(x-xmax)*(xdim/2)/xmax;
  int yimage = -(y-ymax)*(ydim/2)/ymax;
  
  
  DGtal::Z2i::RealPoint NewPoint(ximage,yimage);
  
  
  return NewPoint;
  
}
*/


bool
CoronaryArteryTree::hasIntersections(Segment<Point2D> S1, Point2D newPoint)
{
  DGtal::Z2i::RealPoint Barycenter = findBarycenter(newPoint, S1.myIndex);
  
  Segment<Point2D> newSegment;
  newSegment.myCoordinate = newPoint;
  
  Segment<Point2D> rootSegment;
  rootSegment.myCoordinate = newPoint;
  Segment<Point2D> newOldSegment;
  newOldSegment.myCoordinate = newPoint;
  for (auto s : myVectSegments)
  {
    if(s!=newSegment && s!=rootSegment && s!= newOldSegment)
    {
      
      if(hasIntersection(newPoint ,Barycenter,s.myCoordinate,myVectSegments[myVectParent[s.myIndex]].myCoordinate) || hasIntersection(s.myCoordinate,myVectSegments[myVectParent[s.myIndex]].myCoordinate,Barycenter,S1.myCoordinate) || hasIntersection(s.myCoordinate,myVectSegments[myVectParent[s.myIndex]].myCoordinate,Barycenter,myVectSegments[myVectParent[S1.myIndex]].myCoordinate) )
      {
        return true;
      }
    }
    
  }
  return false;
}

bool operator==(CoronaryArteryTree::Segment<CoronaryArteryTree::Point2D>  S1, CoronaryArteryTree::Segment<CoronaryArteryTree::Point2D> S2)
{
  if(S1.myCoordinate==S2.myCoordinate)
  {
    return true;
  }
  else
  {
    return false;
  }
}

bool operator!=(CoronaryArteryTree::Segment<CoronaryArteryTree::Point2D>  S1, CoronaryArteryTree::Segment<CoronaryArteryTree::Point2D> S2)
{
  if(S1.myCoordinate!=S2.myCoordinate)
  {
    return true;
  }
  else
  {
    return false;
  }
}
