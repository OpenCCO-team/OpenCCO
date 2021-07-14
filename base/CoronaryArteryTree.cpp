
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
CoronaryArteryTree::updateFlowTerminal(unsigned int segIndex)
{
  //Update HydroResistance
  myVectSegments[segIndex].myHydroResistance = 8.0*my_mu*myVectSegments[segIndex].myLength/M_PI;
  assert(myVectSegments[segIndex].myHydroResistance != 0);
  //std::cout<<"myHydroResistance="<<myVectSegments[segIndex].myHydroResistance<<std::endl;
}

void
CoronaryArteryTree::updateFlowParameters(unsigned int segIndex)
{
  if(myVectSegments[segIndex].myIndex != 0) {
    //Update HydroResistance
    myVectSegments[segIndex].myHydroResistance = 8.0*my_mu*myVectSegments[segIndex].myLength/M_PI;
    if(myVectChildren[segIndex].first!=0 && myVectChildren[segIndex].second!=0) {
      // Get left and right children of the segment
      Segment<Point2D> sLeft = myVectSegments[myVectChildren[segIndex].first];
      Segment<Point2D> sRight = myVectSegments[myVectChildren[segIndex].second];
      double r1 = (sLeft.myRadius/myVectSegments[segIndex].myRadius);
      double r2 = (sRight.myRadius/myVectSegments[segIndex].myRadius);
      myVectSegments[segIndex].myHydroResistance += 1.0/((r1*r1*r1*r1)/sLeft.myHydroResistance + (r2*r2*r2*r2)/sRight.myHydroResistance) ;
    }
    assert(myVectSegments[segIndex].myHydroResistance != 0);
    //std::cout<<"myHydroResistance="<<myVectSegments[segIndex].myHydroResistance<<std::endl;
    //Update beta
    if(myVectSegments[myVectParent[segIndex]].myIndex!=0) { //it's not the first segment
      std::pair<int, int> children = myVectChildren[myVectParent[segIndex]];
      int brotherIndex;
      if(children.first != segIndex)
        brotherIndex = children.first;
      else
        brotherIndex = children.second;
      double segFlow = myVectSegments[segIndex].myFlow;
      double brotherFlow = myVectSegments[brotherIndex].myFlow;
      double segHydro = myVectSegments[segIndex].myHydroResistance;
      double brotherHydro = myVectSegments[brotherIndex].myHydroResistance;
      double ratioR = pow((segFlow*segHydro) / (brotherFlow*brotherHydro), 0.25); //ratio of radii of two brother segments
      //ratio of radii of current segment wrt the parent segment
      myVectSegments[segIndex].beta = pow(1.0/(1+pow(ratioR,-my_gamma)),1.0/my_gamma);
      myVectSegments[brotherIndex].beta = pow(1.0/(1+pow(1.0/ratioR,-my_gamma)),1.0/my_gamma);
    }
    updateFlowParameters(myVectSegments[myVectParent[segIndex]].myIndex);
  }
}

void
CoronaryArteryTree::updateRadius(unsigned int segIndex)
{
  if(myVectChildren[segIndex].first!=0 && myVectChildren[segIndex].second!=0) {
    //update radii left and right children
    Segment<Point2D> sLeft = myVectSegments[myVectChildren[segIndex].first];
    Segment<Point2D> sRight = myVectSegments[myVectChildren[segIndex].second];
    double rL = myVectSegments[segIndex].myRadius * sLeft.beta;
    myVectSegments[myVectChildren[segIndex].first].myRadius = rL;
    double rR = myVectSegments[segIndex].myRadius * sRight.beta;
    myVectSegments[myVectChildren[segIndex].second].myRadius = rR;
    updateRadius(sLeft.myIndex);
    updateRadius(sRight.myIndex);
  }
}


void
CoronaryArteryTree::updateRootRadius()
{
  unsigned int idRoot = 1;
  double Qpk = myVectSegments[idRoot].myKTerm*my_qTerm;
  double rRoot = pow((myVectSegments[idRoot].myHydroResistance*Qpk)/(my_pPerf-my_pTerm),0.25);
  myVectSegments[idRoot].myRadius = rRoot;
  updateRadius(idRoot);
}

void
CoronaryArteryTree::updateSegmentRadiusToRoot(unsigned int segIndex)
{
  unsigned int idRoot = 1;
  double rRoot= myVectSegments[idRoot].myRadius;
  std::vector<unsigned int> v = getPathToRoot(myVectSegments[segIndex]);
  double radius = 1;
  for(size_t id=v.size()-1; id>=2; id--) {
    radius *= rRoot*myVectSegments[v[id]].beta;
    myVectSegments[v[id]].myRadius = radius;
  }
}

double
CoronaryArteryTree::computeTreeVolume(double mu, double lambda)
{
  double volume = 0;
  for(size_t idSeg = 1; idSeg<myVectSegments.size(); idSeg++) {
    volume += pow(myVectSegments[idSeg].myLength,mu)*pow(myVectSegments[idSeg].myRadius,lambda);
  }
  return volume;
}

void
CoronaryArteryTree::updateScale(double scale)
{
  for(size_t it=0; it<this->myVectSegments.size(); it++) {
    this->myVectSegments[it].myCoordinate[0] *= scale;
    this->myVectSegments[it].myCoordinate[1] *= scale;
    if(it!=0) //except root radius
      this->myVectSegments[it].myRadius *= scale;
    this->myVectSegments[it].myLength *= scale;
  }
  this->myRsupp *= scale;
  this->my_rPerf = this->myRsupp;
  this->myDThresold=sqrt((M_PI*this->myRsupp*this->myRsupp)/this->myKTerm);
}

double
CoronaryArteryTree::computeTotalVolume(unsigned int segIndex)
{
  //Volumn of current segment
  double v = myVectSegments[segIndex].myRadius*myVectSegments[segIndex].myLength;
  if(myVectChildren[segIndex].first!=0 && myVectChildren[segIndex].second!=0) {
    double vL = computeTotalVolume(myVectChildren[segIndex].first);
    double vR = computeTotalVolume(myVectChildren[segIndex].second);
    v += vL + vR;
  }
  return v;
}

bool
CoronaryArteryTree::isAddable(const Point2D &p, unsigned int segIndex, unsigned int nbIter, unsigned int nbNeibour)
{
  DGtal::Z2i::RealPoint pOpt;
  double r0 = 1.0, r1 = 1.0, r2 = 1.0;
  
  Point2D newCenter = FindBarycenter(p, segIndex);
  // Creation of the left child
  Segment<Point2D> sNewLeft;
  sNewLeft.myCoordinate = myVectSegments[segIndex].myCoordinate;
  sNewLeft.myRadius = myVectSegments[segIndex].myRadius;
  sNewLeft.myIndex = myVectSegments.size();
  sNewLeft.myLength = (newCenter - myVectSegments[segIndex].myCoordinate).norm();
  sNewLeft.myFlow = myVectSegments[segIndex].myFlow;
  
  // Creation of the right child
  Segment<Point2D> sNewRight;
  sNewRight.myCoordinate = p;
  sNewRight.myRadius = myVectSegments[segIndex].myRadius;
  sNewRight.myIndex = myVectSegments.size()+1;
  sNewRight.myLength = (newCenter - p).norm();
  sNewRight.myFlow = my_qTerm;
  
  
  // Cretation a copy of center segment
  Segment<Point2D> sCurrent = myVectSegments[segIndex];
  sCurrent.myCoordinate = newCenter;
  sCurrent.myLength = (newCenter - myVectSegments[myVectParent[segIndex]].myCoordinate).norm();
  //sCurrent.myRadius = myVectSegments[segIndex].myRadius;
  sCurrent.myFlow = myVectSegments[segIndex].myFlow + my_qTerm;
  
  // Cretation a copy of parent segment
  Segment<Point2D> sParent = myVectSegments[myVectParent[segIndex]];
  
  bool res1 = kamyiaOptimization(sParent.myCoordinate, sCurrent, sNewLeft, sNewRight, nbIter, pOpt, r0, r1, r2);
  bool res2;
  if(res1) {
    res2 = isIntersecting(p, pOpt, segIndex, nbNeibour);
    if(!res2) {
      //Update optimal values
      sNewLeft.myLength = (sNewLeft.myCoordinate - pOpt).norm();
      sNewLeft.myRadius = r1;
      sNewRight.myLength = (sNewRight.myCoordinate - pOpt).norm();
      sNewRight.myRadius = r2;
      //sCurrent.myCoordinate = pOpt;
      //sCurrent.myRadius = r0;
      //sCurrent.myLength = (sParent.myCoordinate - pOpt).norm();
      
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
      myVectSegments[segIndex].myRadius = r0;
      myVectSegments[segIndex].myLength = (pOpt - myVectSegments[myVectParent[segIndex]].myCoordinate).norm();
      //update childrens of center segment
      myVectChildren[segIndex].first = sNewLeft.myIndex;
      myVectChildren[segIndex].second = sNewRight.myIndex;
      
      //Update center coordinate
      myVectSegments[segIndex].myCoordinate = pOpt;
      
      // update parameters
      myKTerm++;
      
      // Update physilogique paramaters
      updateFlowTerminal(sNewRight.myIndex);
      updateFlowParameters(sNewLeft.myIndex);
      
      // Update root radius
      updateRootRadius();
      
    }
  }
  return res1 && !res2;
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
  Point2D newCenter = FindBarycenter(p, nearIndex);
  
  // Creation of the left child
  Segment<Point2D> sNewLeft;
  sNewLeft.myCoordinate = myVectSegments[nearIndex].myCoordinate;
  sNewLeft.myRadius = rLeft;
  sNewLeft.myIndex = myVectSegments.size();
  sNewLeft.myLength = (newCenter - myVectSegments[nearIndex].myCoordinate).norm();
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
  sNewRight.myLength = (newCenter - p).norm();
  sNewRight.myFlow = my_qTerm;
  myVectSegments.push_back(sNewRight);
  myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0));
  myVectParent.push_back(nearIndex);
  myVectTerminals.push_back(sNewRight.myIndex);
  
  // Update center segment
  myVectSegments[nearIndex].myCoordinate = newCenter;
  myVectSegments[nearIndex].myLength = (newCenter - myVectSegments[myVectParent[nearIndex]].myCoordinate).norm();
  //update childrens of center segment
  myVectChildren[nearIndex].first = sNewLeft.myIndex;
  myVectChildren[nearIndex].second = sNewRight.myIndex;
  /*
  //Avant
  std::string filename = "testCCO_"+std::to_string(nearIndex)+"A.eps";
  exportBoardDisplay(filename.c_str(), true);
  myBoard.clear();
  */
  // Optimization
  bool res = kamyiaOptimization(nearIndex);
  /*
  //Apres
  filename = "testCCO_"+std::to_string(nearIndex)+"B.eps";
  exportBoardDisplay(filename.c_str(), true);
  myBoard.clear();
  */
  // Update Kterm
  myKTerm++;
  
  // Update physilogique paramaters
  updateFlowTerminal(sNewRight.myIndex);
  updateFlowParameters(sNewLeft.myIndex);
  
  // Update root radius
  updateRootRadius();
  
  //TODO: then opt with volume of the
  
  return true;
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

bool
CoronaryArteryTree::addSegmentFromPointBK(const Point2D &p,
                                        unsigned int nearIndex,
                                        double rLeft, double rRight)
{

  // To add a new segment we need to add two new points (p and middle point) with two segments.
  // (a) First point added : point associated to the nearest segment. (basic solution the center of the nearest)
  // (b) Second point added: the new point with its new segment.
  // to process (a): s
  Point2D newCenter = FindBarycenter(p, nearIndex);
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
  sNewLeft.myLength = (newCenter - myVectSegments[nearIndex].myCoordinate).norm();
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
  sNewRight.myLength = (newCenter - p).norm();
  myVectSegments.push_back(sNewRight);
  myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0));
  myVectParent.push_back(nearIndex);
  myVectTerminals.push_back(sNewRight.myIndex);
  // Update center segment
  myVectSegments[nearIndex].myCoordinate = newCenter;
  myVectSegments[nearIndex].myLength = (newCenter - myVectSegments[myVectParent[nearIndex]].myCoordinate).norm();
  //update childrens of center segment
  myVectChildren[nearIndex].first = sNewLeft.myIndex;
  myVectChildren[nearIndex].second = sNewRight.myIndex;

  myKTerm++;

  return true;
}



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
                                       bool updateDisplay ){
  if (updateDisplay){
    boardDisplay(thickness);
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


CoronaryArteryTree::Point2D
CoronaryArteryTree::getSegmentCenter(unsigned int i)
{
  return getSegmentCenter(myVectSegments.at(i));
}


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



unsigned int
CoronaryArteryTree::getParentSegment(const Segment<Point2D> &s){
  return myVectParent[s.myIndex];
}

unsigned int
CoronaryArteryTree::getLeftChild(const Segment<Point2D> &s){
  return myVectChildren[s.myIndex].first;
}

unsigned int
CoronaryArteryTree::getRightChild(const Segment<Point2D> &s){
  return myVectChildren[s.myIndex].second;
}

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


void
CoronaryArteryTree::udpatePerfusionArea(){
  myCurrAPerf = (myKTerm+1.0)*M_PI*myRsupp*myRsupp;
  
}

CoronaryArteryTree::Point2D
CoronaryArteryTree::generateNewLocation(unsigned int nbTrials){
  Point2D res;
  bool found = false;
  unsigned int n = nbTrials;
  while(!found && n >= 0) {
    n--;
    auto p = generateALocation();
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
CoronaryArteryTree::generateALocation() {
  Point2D res;
  res = generateRandomPtOnDisk(myTreeCenter, myRsupp);
  bool isComp = true;
  unsigned int id = 1;
  while ( isComp && id < myVectTerminals.size() ) {
    
    isComp = (myVectSegments[myVectTerminals[id]].myCoordinate - res).norm() > myDThresold;
    id++;
  }
  return  std::pair<Point2D, bool> {res, isComp};
}

double
CoronaryArteryTree::getDistance(unsigned int index,
                                const Point2D &p ) const {
  return (myVectSegments[index].myCoordinate-p).norm();
}

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
  double ratioQ = sL.myFlow/(myVectSegments[index].myFlow+my_qTerm);//0.5;
  std::cout<<"ratioQ="<< sL.myFlow/myVectSegments[index].myFlow<<std::endl;
  double r0 = myVectSegments[index].myRadius;
  double R0 = r0*r0; //R0 = r0*r0
  double R1 = r0*r0; //R1 = r1*r1
  double R2 = r0*r0; //R2 = r2*r2
  
  double f0 =  R0*R0*R0;
  double f1 = ratioQ*f0;//k * r1*r1*r1; //left
  double f2 = (1.0-ratioQ)*f0;//k * r2*r2*r2; //right
  // Starting position from Equation (21) for initialisation as mentionned page 11 [Clara Jaquet et HT]
  DGtal::Z2i::RealPoint pb ((f0*pParent[0]+f1*pL[0]+f2*pR[0])/(2.0*f0), (f0*pParent[1]+f1*pL[1]+f2*pR[1])/(2.0*f0));
  std::cout<<"Init:"<<pb<<std::endl;
  double l0 = (pParent - pb).norm();
  double l1 = (pL - pb).norm();
  double l2 = (pR - pb).norm();
  //std::cout<<"l0="<<l0<<", l1="<<l1<<", l2="<<l2<<std::endl;
  double deltaP1 = (f0*l0)/(R0*R0)+(f1*l1)/(R1*R1);
  double deltaP2 = (f0*l0)/(R0*R0)+(f2*l2)/(R2*R2);
  
  double rr1 = sParent.myRadius;
  double rr2 = sParent.myRadius;
  // The variable to solve for with its initial value. It will be
  // mutated in place by the solver.
  //DGtal::Z2i::RealPoint pbInit ((f0*pParent[0]+f1*pL[0]+f2*pR[0])/(2.0*f0), (f0*pParent[1]+f1*pL[1]+f2*pR[1])/(2.0*f0));
  //DGtal::trace.info() << pbInit << std::endl;
  
  bool hasSolution = true;
  
  for (int i = 0; i< nbIter && hasSolution; i++){
    DGtal::trace.progressBar(i, nbIter);
    l0 = (pParent - pb).norm();
    l1 = (pL - pb).norm();
    l2 = (pR - pb).norm();
    
    hasSolution = kamiyaOpt(deltaP1, deltaP2, f0, f1, f2, l0, l1, l2, rr1, rr2);
    // Equation 27
    R0 = pow(f0*(pow(rr1, my_gamma)/f1 + pow(rr2, my_gamma)/f2), 1.0/my_gamma);
    // Equation (26) page 13
    pb[0] = (pParent[0]*R0/l0 + pL[0]*R1/l1 + pR[0]*R2/l2)/(R0/l0+R1/l1+R2/l2);
    pb[1] = (pParent[1]*R0/l0 + pL[1]*R1/l1 + pR[1]*R2/l2)/(R0/l0+R1/l1+R2/l2);
    deltaP1 = (f0*l0)/(R0*R0)+(f1*l1)/(R1*R1);
    deltaP2 = (f0*l0)/(R0*R0)+(f2*l2)/(R2*R2);
    std::cout << "xpL[0] : " << rr1 << " and xpR[0] " << rr2 << " r0 " << R0 << "\n";
    std::cout << "pbNew[0]  : " << pb[0] << " and pbNew[1] " << pb[1] << "\n";
  }
  // If there is a solution, then update result
  if(hasSolution) {
    // Update center point coordinate
    //myVectSegments[myVectParent[index]].myCoordinate[0] = pb[0];
    //myVectSegments[myVectParent[index]].myCoordinate[1] = pb[1];
    myVectSegments[index].myCoordinate[0] = pb[0];
    myVectSegments[index].myCoordinate[1] = pb[1];
    double mR0 = sqrt(R0);
    double mR1 = sqrt(rr1);
    double mR2 = sqrt(rr2);
    double mL0 = (pb - pParent).norm();
    double mL1 = (pb - pL).norm();
    double mL2 = (pb - pR).norm();
    
    // Update radius of center, left and right segments
    myVectSegments[index].myRadius = mR0;
    myVectSegments[myVectChildren[index].first].myRadius =  mR1;
    myVectSegments[myVectChildren[index].second].myRadius =  mR2;
    // Update length segments
    myVectSegments[index].myLength = mL0;
    myVectSegments[myVectChildren[index].first].myLength = mL1;
    myVectSegments[myVectChildren[index].second].myLength = mL2;
    
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
  //Segment<Point2D> sParent = myVectSegments[myVectParent[index]];
  //Segment<Point2D> sL = myVectSegments[myVectChildren[index].first];
  //Segment<Point2D> sR = myVectSegments[myVectChildren[index].second];
  
  //DGtal::Z2i::RealPoint pParent = sParent.myCoordinate;
  DGtal::Z2i::RealPoint pL = sL.myCoordinate;
  DGtal::Z2i::RealPoint pR = sR.myCoordinate;
  
  //std::cout<<"myFlow="<< myVectSegments[index].myFlow<<" and sL.myFlow="<<sL.myFlow<<" and sR.myFlow="<<sR.myFlow<<std::endl;
  double ratioQ = sL.myFlow/sCurrent.myFlow;//0.5;
  //std::cout<<"ratioQ="<< sL.myFlow/sCurrent.myFlow<<std::endl;
  double rr0 = sCurrent.myRadius;
  double R0 = rr0*rr0; //R0 = r0*r0
  double R1 = rr0*rr0; //R1 = r1*r1
  double R2 = rr0*rr0; //R2 = r2*r2
  
  double f0 =  rr0*rr0*rr0;
  double f1 = ratioQ*f0;//k * r1*r1*r1; //left
  double f2 = (1.0-ratioQ)*f0;//k * r2*r2*r2; //right
  // Starting position from Equation (21) for initialisation as mentionned page 11 [Clara Jaquet et HT]
  DGtal::Z2i::RealPoint pb ((f0*pParent[0]+f1*pL[0]+f2*pR[0])/(2.0*f0), (f0*pParent[1]+f1*pL[1]+f2*pR[1])/(2.0*f0));
  //std::cout<<"Init:"<<pb<<std::endl;
  double l0 = (pParent - pb).norm();
  double l1 = (pL - pb).norm();
  double l2 = (pR - pb).norm();
  //std::cout<<"l0="<<l0<<", l1="<<l1<<", l2="<<l2<<std::endl;
  double deltaP1 = (f0*l0)/(R0*R0)+(f1*l1)/(R1*R1);
  double deltaP2 = (f0*l0)/(R0*R0)+(f2*l2)/(R2*R2);
  
  double rr1 = sCurrent.myRadius;
  double rr2 = sCurrent.myRadius;
  // The variable to solve for with its initial value. It will be
  // mutated in place by the solver.
  //DGtal::Z2i::RealPoint pbInit ((f0*pParent[0]+f1*pL[0]+f2*pR[0])/(2.0*f0), (f0*pParent[1]+f1*pL[1]+f2*pR[1])/(2.0*f0));
  //DGtal::trace.info() << pbInit << std::endl;
  
  bool hasSolution = true;
  
  for (int i = 0; i< nbIter && hasSolution; i++){
    //DGtal::trace.progressBar(i, nbIter);
    l0 = (pParent - pb).norm();
    l1 = (pL - pb).norm();
    l2 = (pR - pb).norm();
    
    hasSolution = kamiyaOpt(deltaP1, deltaP2, f0, f1, f2, l0, l1, l2, R1, R2);
    // Equation 27
    R0 = pow(f0*(pow(R1, my_gamma)/f1 + pow(R2, my_gamma)/f2), 1.0/my_gamma);
    // Equation (26) page 13
    pb[0] = (pParent[0]*R0/l0 + pL[0]*R1/l1 + pR[0]*R2/l2)/(R0/l0+R1/l1+R2/l2);
    pb[1] = (pParent[1]*R0/l0 + pL[1]*R1/l1 + pR[1]*R2/l2)/(R0/l0+R1/l1+R2/l2);
    deltaP1 = (f0*l0)/(R0*R0)+(f1*l1)/(R1*R1);
    deltaP2 = (f0*l0)/(R0*R0)+(f2*l2)/(R2*R2);
    //std::cout << "xpL[0] : " << rr1 << " and xpR[0] " << rr2 << " r0 " << R0 << "\n";
    //std::cout << "pbNew[0]  : " << pb[0] << " and pbNew[1] " << pb[1] << "\n";
  }
  // If there is a solution, then update result
  if(hasSolution) {
    pOpt = pb;
    r0 = sqrt(R0);
    r1 = sqrt(R1);
    r2 = sqrt(R2);
    double mL0 = (pb - pParent).norm();
    double mL1 = (pb - pL).norm();
    double mL2 = (pb - pR).norm();
    
    if(2*r0<=mL0 && 2*r1<=mL1 && 2*r2<mL2)
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
  out << "main parameters: myKTerm: "  << myKTerm << "\n\t myDThresold: " << myDThresold << "\n\t myRsupp: " << myRsupp << std::endl;
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
  Point2D newCenter = FindBarycenter(p, nearIndex);
  
  // Creation of the left child
  Segment<Point2D> sNewLeft;
  sNewLeft.myCoordinate = myVectSegments[nearIndex].myCoordinate;
  sNewLeft.myRadius = 1.0;
  sNewLeft.myIndex = myVectSegments.size();
  sNewLeft.myLength = (newCenter - myVectSegments[nearIndex].myCoordinate).norm();
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
  sNewRight.myLength = (newCenter - p).norm();
  sNewRight.myFlow = my_qTerm;
  myVectSegments.push_back(sNewRight);
  myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0,0));
  myVectParent.push_back(nearIndex);
  myVectTerminals.push_back(sNewRight.myIndex);
  
  // Update center segment
  myVectSegments[nearIndex].myCoordinate = newCenter;
  myVectSegments[nearIndex].myLength = (newCenter - myVectSegments[myVectParent[nearIndex]].myCoordinate).norm();
  //update childrens of center segment
  myVectChildren[nearIndex].first = sNewLeft.myIndex;
  myVectChildren[nearIndex].second = sNewRight.myIndex;
  
  // Update Kterm
  myKTerm++;
  
  kamyiaOptimization(nearIndex);
  
  //updateRadius();
  
  return true;
  
}

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

void
CoronaryArteryTree::updateTreshold()
{
  myDThresold=sqrt((M_PI*myRsupp*myRsupp)/myKTerm);
}
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







double
CoronaryArteryTree::GetLength(unsigned int Index)
{
  Point2D distalPoint=myVectSegments[Index].myCoordinate;
  Point2D proximalPoint = myVectSegments[myVectParent[Index]].myCoordinate;
  double length  = pow(pow((proximalPoint[0] - distalPoint[0]), 2.0) + pow((proximalPoint[1] - distalPoint[1]), 2.0),0.5);
  return length;
}



double
CoronaryArteryTree::GetTotalVolume(const Point2D &p1,const Point2D &p2,const Point2D &p3,const Point2D &pOpti)
{
  double H1,H2,H3;
  
  H1 = sqrt(pow(p1[0]-pOpti[0],2)+pow(p1[1]-pOpti[1],2));
  H2 = sqrt(pow(p2[0]-pOpti[0],2)+pow(p2[1]-pOpti[1],2));
  H3 = sqrt(pow(p3[0]-pOpti[0],2)+pow(p3[1]-pOpti[1],2));
  
  return M_PI*H1+M_PI*H2+M_PI*H3;
  
}

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
CoronaryArteryTree::FindBarycenter(const Point2D &p, unsigned int index)
{
  Point2D first_point=p;
  Point2D second_point=myVectSegments[index].myCoordinate;
  Point2D third_point=myVectSegments[myVectParent[index]].myCoordinate;
  Point2D barycenter((first_point[0]+second_point[0]+third_point[0])/3.0,(first_point[1]+second_point[1]+third_point[1])/3.0);
  return barycenter;
}


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

double
CoronaryArteryTree::FindXmax(int xDim, int yDim)
{
  
  return myRsupp*cos(atan(xDim/yDim));
}

double
CoronaryArteryTree::FindYmax(int xDim, int yDim)
{
  
  return myRsupp*sin(atan(xDim/yDim));
  
}


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



bool
CoronaryArteryTree::hasIntersections(Segment<Point2D> S1, Point2D newPoint)
{
  DGtal::Z2i::RealPoint Barycenter = FindBarycenter(newPoint, S1.myIndex);
  
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
