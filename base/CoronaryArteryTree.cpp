
#include "CoronaryArteryTree.h"
#include "geomhelpers.h"
#include "ConstructionHelpers.h"

#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/readers/PPMReader.h"
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/colormaps/GradientColorMap.h"

#include "DGtal/images/ArrayImageIterator.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/geometry/helpers/ContourHelper.h"

#include <math.h>

CoronaryArteryTree::Segment<CoronaryArteryTree::Point2D>
CoronaryArteryTree::updateResistanceFromRoot(unsigned int segIndex) {
  if(segIndex != 0) {
    //Update HydroResistance
    double l = getLengthSegment(segIndex);
    myVectSegments[segIndex].myResistance = 8.0*my_mu*l/M_PI;
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

void
CoronaryArteryTree::updateFlow(unsigned int segIndex)
{
  Segment<Point2D> sC = myVectSegments[segIndex];
  while (sC.myIndex != 0){
    myVectSegments[sC.myIndex].myKTerm ++;
    myVectSegments[sC.myIndex].myFlow = myVectSegments[sC.myIndex].myKTerm*my_qTerm;
    sC = myVectSegments[myVectParent[sC.myIndex]];
  }
}

void
CoronaryArteryTree::updateLengthFactor()
{
  double r_pk = sqrt(myKTerm + 1)*myRsupp;
  myLengthFactor =  r_pk / my_rPerf;
}

void
CoronaryArteryTree::updateRadius(unsigned int segIndex, double beta)
{
  if(myVectChildren[segIndex].first!=0 && myVectChildren[segIndex].second!=0) {
    //update radii left and right children
    Segment<Point2D> sLeft = myVectSegments[myVectChildren[segIndex].first];
    Segment<Point2D> sRight = myVectSegments[myVectChildren[segIndex].second];
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
  double Qpk = myVectSegments[idRoot].myKTerm*my_qTerm;
  double rRoot = pow((myVectSegments[idRoot].myResistance*Qpk)/my_pDrop,0.25);
  myVectSegments[idRoot].myRadius = rRoot;
  updateRadius(idRoot, 1.0);
}

double
CoronaryArteryTree::getLengthSegment(unsigned int segIndex)
{
  if(segIndex==0) return 0; //there is no parent segment
  Point2D p1 = myVectSegments[segIndex].myCoordinate;
  Point2D p2 = myVectSegments[myVectParent[segIndex]].myCoordinate;
  return (p1-p2).norm()*myLengthFactor;
}

double
CoronaryArteryTree::computeTotalVolume(unsigned int segIndex)
{
  //Volumn of current segment : Eq 20
  double l = getLengthSegment(segIndex);
  double r = myVectSegments[segIndex].myRadius;
  double v = myVectSegments[segIndex].myRadius*myVectSegments[segIndex].myRadius*getLengthSegment(segIndex);
  double vL, vR;
  if(myVectChildren[segIndex].first!=0 && myVectChildren[segIndex].second!=0) {
    vL = computeTotalVolume(myVectChildren[segIndex].first);
    vR = computeTotalVolume(myVectChildren[segIndex].second);
    v += vL + vR;
  }
  return v;
}

bool
CoronaryArteryTree::isAddable(const Point2D &p, unsigned int segIndex,
                              unsigned int nbIter, double tolerance,
                              unsigned int nbNeibour, bool verbose)
{
  Point2D pOpt, newCenter = findBarycenter(p, segIndex);
  if(myIsImageDomainRestrained && (myImageDomain(DGtal::Z2i::Point(static_cast<int>(newCenter[0]),
                                      static_cast<int>(newCenter[1])))< 128)){
    return false;
  }
  if(myIsImageDomainRestrained && (!checkNoIntersectDomain(myImageDomain, 128,
                  DGtal::Z2i::Point(static_cast<int>(p[0]), static_cast<int>(p[1])),
                  DGtal::Z2i::Point(static_cast<int>(newCenter[0]),static_cast<int>(newCenter[1])))))
     {
      return false;
  }
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
  Point2D pParent = sParent.myCoordinate;
  /*
  double f1 = sNewLeft.myFlow;//ratioQ*f0;//k * r1*r1*r1; //left
  double f2 = sNewRight.myFlow;//(1.0-ratioQ)*f0;//k * r2*r2*r2; //right
  double f0 =  f1 + f2;//R0*R0*R0; //middle
  Point2D pL = sNewLeft.myCoordinate;
  Point2D pR = sNewRight.myCoordinate;
  Point2D pCurrent ((f0*pParent[0]+f1*pL[0]+f2*pR[0])/(2.0*f0), (f0*pParent[1]+f1*pL[1]+f2*pR[1])/(2.0*f0)); //barycenter of the sement
  */
  Point2D pCurrent = (sParent.myCoordinate+sNewLeft.myCoordinate)/2.0; //center of the current segment
  double r0 = sCurrent.myRadius, r1 = sCurrent.myRadius, r2 = sCurrent.myRadius;
  bool res1 = true, res2 = false, res3 = true, isDone = false;
  double vol = -1, volCurr = -1, diffVol = -1;
  CoronaryArteryTree cTreeCurr = *this;
  if (verbose){
    DGtal::trace.info() <<"---------- segIndex: "<<segIndex<<std::endl;
  }
  size_t i=0;
  for(size_t i=0; i<nbIter && !isDone; i++) {
    res1 = kamyiaOptimization(pCurrent, pParent, sCurrent.myRadius, sNewLeft, sNewRight, 1, pOpt, r0, r1, r2);
    if(!res1) //Kamyia does not have solution
      return false;
    
    res2 = isIntersecting(p, pOpt, segIndex, nbNeibour, 2*cTreeCurr.myVectSegments[segIndex].myRadius);
    if(res2) //there is intersection betweeen segments
      return false;
    
    if(myIsImageDomainRestrained && (myImageDomain(DGtal::Z2i::Point(static_cast<int>(pOpt[0]), static_cast<int>(pOpt[1])))< 128)) //resulting point is out of the domaine
      return false;
    
    //resulting segments intersect the domaine
    if(myIsImageDomainRestrained && (!checkNoIntersectDomain(myImageDomain, 128,
                    DGtal::Z2i::Point(static_cast<int>(p[0]), static_cast<int>(p[1])),
                    DGtal::Z2i::Point(static_cast<int>(pOpt[0]),static_cast<int>(pOpt[1]))))) //Center segment
      return false;
    if(myIsImageDomainRestrained && (!checkNoIntersectDomain(myImageDomain, 128,
                    DGtal::Z2i::Point(static_cast<int>(sNewLeft.myCoordinate[0]), static_cast<int>(sNewLeft.myCoordinate[1])),
                    DGtal::Z2i::Point(static_cast<int>(pOpt[0]),static_cast<int>(pOpt[1]))))) //Left segment
      return false;
    if(myIsImageDomainRestrained && (!checkNoIntersectDomain(myImageDomain, 128,
                    DGtal::Z2i::Point(static_cast<int>(sNewRight.myCoordinate[0]), static_cast<int>(sNewRight.myCoordinate[1])),
                    DGtal::Z2i::Point(static_cast<int>(pOpt[0]),static_cast<int>(pOpt[1]))))) //Right segment
      return false;
    
    //Otherwise, iterate for the optimisation process
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
    if (verbose){
      DGtal::trace.info() <<"Iter "<<i<<" has tree Volume: "<< vol <<std::endl;
    }
    if(i==0) {
      volCurr = vol;
      cTreeCurr = cTree1;
    }
    else {
      diffVol = volCurr - vol;
      if(fabs(diffVol) < tolerance) {
        // Verify the degenerate case of the resulting segment
        double l0 = (pOpt - pParent).norm()*myLengthFactor;
        double l1 = (pOpt - sNewLeft.myCoordinate).norm()*myLengthFactor;
        double l2 = (pOpt - sNewRight.myCoordinate).norm()*myLengthFactor;
        res3 = (2*r0<=l0) && (2*r1<=l1) && (2*r2<=l2);
        // If there is a solution, then save the result
        if(res3) {
          //std::cout<<"Best volume at "<<i<<std::endl;
          isDone = true;
          *this = cTree1;
        }
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
  
  return isDone;
}


bool
CoronaryArteryTree::addSegmentFromPoint(const Point2D &p,
                                        unsigned int nearIndex)
{

  Point2D newCenter = findBarycenter(p, nearIndex);
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
  sNewLeft.myIndex = myVectSegments.size();
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
  Segment<Point2D> sNewRight;
  sNewRight.myCoordinate = p;
  sNewRight.myRadius = myVectSegments[nearIndex].myRadius;
  sNewRight.myIndex = myVectSegments.size();
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

void
CoronaryArteryTree::boardDisplay(double thickness, bool clearDisplay)
{
  double scaleFactorEPS = 0.01;
  myBoard.setUnit(scaleFactorEPS, LibBoard::Board::UCentimeter);
  if (clearDisplay){
    myBoard.clear();
  }
  // 57.5 from myBoard change scale
  double scaleBoard = 57.5*scaleFactorEPS;
  // drawing base circle
  if (!myIsImageDomainRestrained){
    myBoard.setPenColor(DGtal::Color::Blue);
    myBoard.setLineWidth(myVectSegments[0].myRadius*scaleBoard*thickness);
    myBoard.drawCircle(myTreeCenter[0], myTreeCenter[1], my_rPerf, 1);
  }
  
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
      continue;
  
    // distal node
    Point2D distal = s.myCoordinate;
    Point2D proxital = myVectSegments[myVectParent[s.myIndex]].myCoordinate;
    myBoard.setLineWidth(myVectSegments[s.myIndex].myRadius*scaleBoard*thickness);
    // myBoard.setPenColor(cmap_grad(i));
    myBoard.setPenColor(DGtal::Color::Gray);
    myBoard.drawLine(distal[0], distal[1], proxital[0], proxital[1],2);
    myBoard.setLineWidth(1.0);
    myBoard.setPenColor(DGtal::Color(180, 0, 0, 180));
    myBoard.fillCircle(distal[0], distal[1], myVectSegments[s.myIndex].myRadius*thickness, 1 );
    i++;
  }
  if (myIsImageDomainRestrained){
    std::vector<std::vector<DGtal::Z2i::Point>> vectContours = ConstructionHelpers::getImageContours(myImageDomain, 128);
    for(auto const c: vectContours){
      DGtal::Color col;
      if(DGtal::ContourHelper::isCounterClockWise(c)){
        col = DGtal::Color(200, 200, 200);
      }else{
        col = DGtal::Color::White;
      }
      std::vector<LibBoard::Point> bv;
      for(unsigned int i=0; i<c.size(); i++){
        bv.push_back(LibBoard::Point(c[i][0], c[i][1]));
      }
      myBoard.setFillColor(col);
      myBoard.setLineWidth(5.0*scaleFactorEPS);
      myBoard.fillPolyline(bv);
      myBoard.drawPolyline(bv);
    }
    
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
  Point2D res;
  if (myIsImageDomainRestrained){
    res = generateRandomPtOnImageDomain<CoronaryArteryTree::Point2D>(myImageDomain, 128);
  } else {
    res = generateRandomPtOnDisk(myTreeCenter, my_rPerf);
  }
  bool isComp = true;
  unsigned int id = 1;
  /*
  while ( isComp && id < myVectTerminals.size() ) {
    //isComp = (myVectSegments[myVectTerminals[id]].myCoordinate - res).norm() > myDThresold;
   isComp = getProjDistance(myVectTerminals[id], res) > myDThresold;
    id++;
  }
  */
  //generated point must have a certain distance to ALL tree segments
  while ( isComp && id < myVectSegments.size() ) {
    isComp = getProjDistance(myVectSegments[id].myIndex, res) > myDThresold;
    id++;
  }
  return  std::pair<Point2D, bool> {res, isComp};
}

bool
CoronaryArteryTree::isToCloseFromNearest(const Point2D &p, double minDist) const{
  double d = getProjDistance(getN_NearestSegments(p,1)[0],p);
  return d < minDist;
}
double
CoronaryArteryTree::getProjDistance(const Point2D &p0, const Point2D &p1, const Point2D &p) const{
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
    if (myVectParent[s] == indexPFather && s == indexPChild)
      continue;

    bool inter1 = hasIntersection(p0, pBifurcation, myVectSegments[s].myCoordinate,
                                  myVectSegments[myVectParent[s]].myCoordinate  );
    bool inter2 = hasIntersection(p1, pBifurcation, myVectSegments[s].myCoordinate,
                                  myVectSegments[myVectParent[s]].myCoordinate  );
    bool inter3 = hasIntersection(pAdded, pBifurcation, myVectSegments[s].myCoordinate,
                                  myVectSegments[myVectParent[s]].myCoordinate  );

    if (inter1 || inter2 || inter3)
      return  true;
  }
  return false;
  
}

bool
CoronaryArteryTree::kamyiaOptimization(const Point2D& pCurrent,
                                       const Point2D& pParent,
                                       double rCurrent,
                                       const Segment<Point2D>& sL,
                                       const Segment<Point2D>& sR,
                                       unsigned int nbIter,
                                       Point2D& pOpt,
                                       double& r0, double& r1, double& r2){
  
  // Prepare
  Point2D pL = sL.myCoordinate;
  Point2D pR = sR.myCoordinate;
  
  //std::cout<<"myFlow="<< myVectSegments[index].myFlow<<" and sL.myFlow="<<sL.myFlow<<" and sR.myFlow="<<sR.myFlow<<std::endl;
  //double ratioQ = sL.myFlow/sCurrent.myFlow;//0.5;
  //std::cout<<"ratioQ="<< ratioQ<<std::endl;
  //double rr0 = sCurrent.myRadius;
  double R0 = rCurrent*rCurrent;//rr0*rr0; //R0 = r0*r0
  double R1 = sL.myRadius*sL.myRadius;//rr0*rr0; //R1 = r1*r1
  double R2 = sR.myRadius*sR.myRadius;//rr0*rr0; //R2 = r2*r2
  
  double f1 = sL.myFlow;//ratioQ*f0;//k * r1*r1*r1; //left
  double f2 = sR.myFlow;//(1.0-ratioQ)*f0;//k * r2*r2*r2; //right
  double f0 =  f1 + f2;//R0*R0*R0; //middle
  // Starting position from Equation (21) for initialisation as mentionned page 11 [Clara Jaquet et HT]
  //Point2D pb ((f0*pParent[0]+f1*pL[0]+f2*pR[0])/(2.0*f0), (f0*pParent[1]+f1*pL[1]+f2*pR[1])/(2.0*f0));
  //Point2D pb ((pParent[0]+pL[0])/(2.0), (pParent[1]+pL[1])/(2.0)); //center of the old segment
  Point2D pb (pCurrent[0], pCurrent[1]);
  //std::cout<<"Init:"<<pb<<std::endl;
  double l0 = (pParent - pb).norm()*myLengthFactor;
  double l1 = (pL - pb).norm()*myLengthFactor;
  double l2 = (pR - pb).norm()*myLengthFactor;
  //std::cout<<"l0="<<l0<<", l1="<<l1<<", l2="<<l2<<std::endl;
  double kappa = 8 * my_mu / M_PI;
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
    hasSolution = kamyiaOpt(my_gamma, deltaP1, deltaP2, f0, f1, f2, l0, l1, l2,rr1, rr2);
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

double
CoronaryArteryTree::getDistanceThreshold()
{
  // Check with python code, here is the equation after simplification
  return sqrt((M_PI*my_rPerf*my_rPerf)/myKTerm);
  //Scale factor must be taken into account here
  //return sqrt((M_PI*(myKTerm + 1)*myRsupp*myRsupp)/myKTerm) / myLengthFactor;
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



bool
CoronaryArteryTree::restrainDomain(const std::string &imageName, unsigned int threshold){
  myImageFileDomain = imageName;
  myImageDomain = DGtal::GenericReader<Image>::import( imageName );
  bool isOk = false;
  // Check if at least one pixel of with foreground value exist:
  for (auto p: myImageDomain.domain()){
    if (myImageDomain(p) >= threshold){
      isOk = true;
      break;
    }
  }
  //Check if the root is inside the domaine
  if(isOk) {
    Point2D pRoot = myVectSegments[0].myCoordinate;
    if(myImageDomain(DGtal::Z2i::Point(static_cast<int>(pRoot[0]), static_cast<int>(pRoot[1]))) >= threshold) {
      myIsImageDomainRestrained = true;
      return true;
    }
  }
  myImageFileDomain = "";
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
