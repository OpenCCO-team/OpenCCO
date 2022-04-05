
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

CoronaryArteryTree::Segment<CoronaryArteryTree::Point3D>
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
    Segment<Point3D> sLeft = updateResistanceFromRoot(segIndexLeft);
    Segment<Point3D> sRight = updateResistanceFromRoot(segIndexRight);
    
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
  Segment<Point3D> sC = myVectSegments[segIndex];
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
    Segment<Point3D> sLeft = myVectSegments[myVectChildren[segIndex].first];
    Segment<Point3D> sRight = myVectSegments[myVectChildren[segIndex].second];
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
  Point3D p1 = myVectSegments[segIndex].myCoordinate;
  Point3D p2 = myVectSegments[myVectParent[segIndex]].myCoordinate;
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
CoronaryArteryTree::isAddable(const Point3D &p, unsigned int segIndex, unsigned int nbIter, double tolerance, unsigned int nbNeibour)
{
  Point3D pOpt, newCenter = findBarycenter(p, segIndex);
  // Creation of the left child
  Segment<Point3D> sNewLeft;
  sNewLeft.myCoordinate = myVectSegments[segIndex].myCoordinate;
  sNewLeft.myIndex = myVectSegments.size();
  sNewLeft.myRadius = myVectSegments[segIndex].myRadius;
  sNewLeft.myFlow = myVectSegments[segIndex].myFlow;
  sNewLeft.myKTerm = myVectSegments[segIndex].myKTerm;
  sNewLeft.myResistance = myVectSegments[segIndex].myResistance;
  // Creation of the right child
  Segment<Point3D> sNewRight;
  sNewRight.myCoordinate = p;
  sNewRight.myIndex = myVectSegments.size()+1;
  sNewRight.myRadius = myVectSegments[segIndex].myRadius;
  sNewRight.myFlow = my_qTerm;
  sNewRight.myKTerm = 1;
  // Cretation a copy of center segment
  Segment<Point3D> sCurrent = myVectSegments[segIndex];
  sCurrent.myCoordinate = newCenter;
  sCurrent.myFlow = myVectSegments[segIndex].myFlow + my_qTerm;
  sCurrent.myKTerm = sNewLeft.myKTerm + 1;
  // Cretation a copy of parent segment
  Segment<Point3D> sParent = myVectSegments[myVectParent[segIndex]];
  Point3D pParent = sParent.myCoordinate;
  /*
  double f1 = sNewLeft.myFlow;//ratioQ*f0;//k * r1*r1*r1; //left
  double f2 = sNewRight.myFlow;//(1.0-ratioQ)*f0;//k * r2*r2*r2; //right
  double f0 =  f1 + f2;//R0*R0*R0; //middle
  Point3D pL = sNewLeft.myCoordinate;
  Point3D pR = sNewRight.myCoordinate;
  Point3D pCurrent ((f0*pParent[0]+f1*pL[0]+f2*pR[0])/(2.0*f0), (f0*pParent[1]+f1*pL[1]+f2*pR[1])/(2.0*f0)); //barycenter of the sement
  */
  Point3D pCurrent = (sParent.myCoordinate+sNewLeft.myCoordinate)/2.0; //center of the current segment
  double r0 = sCurrent.myRadius, r1 = sCurrent.myRadius, r2 = sCurrent.myRadius;
  bool res1 = true, res2 = false, res3 = true, isDone = false, res4 = false;
  double vol = -1, volCurr = -1, diffVol = -1;
  CoronaryArteryTree cTreeCurr = *this;
  //std::cout<<"---------- segIndex: "<<segIndex<<std::endl;
  size_t i=0;
  for(size_t i=0; i<nbIter && res1 && !res2 && res3 && !res4 && !isDone; i++) {
    res1 = kamyiaOptimization(pCurrent, pParent, sCurrent.myRadius, sNewLeft, sNewRight, 1, pOpt, r0, r1, r2);
    if(res1) {
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
        
        //Test intersection wrt radius except the childs and the parent and brother
        res4 = cTree1.hasIntersections(sNewLeft.myIndex) || cTree1.hasIntersections(sNewRight.myIndex) || cTree1.hasIntersections(segIndex); //FIXME: epsilon
        if(!res4) {
          vol = cTree1.computeTotalVolume();
          //std::cout<<"Iter "<<i<<" has tree Volume: "<< vol <<std::endl;
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
      }
    }
  }
  return isDone && res1 && !res2 && !res4;
}

bool
CoronaryArteryTree::isIntersecting(const Point3D &pNew, const Point3D &pCenter, unsigned int nearIndex, unsigned int nbNeibour, double minDistance) const
{
  //bool inter = hasNearestIntersections(p, newCenter, 10);
  //bool inter = hasNearestIntersections(myVectParent[nearIndex], nearIndex, p, newCenter,  nbNeibour);
  bool inter = hasIntersections(pNew) || hasIntersections(pCenter);//Old: hasNearestIntersections(myVectParent[nearIndex], nearIndex, pNew, pCenter,  nbNeibour);
  if (inter){
    //DGtal::trace.warning() << "detection intersection" << std::endl;
    return true;
  }
  // Check new point with new segment of barycenter:
  // - newPt and new segment [Barycenter-OriginNewSeg]
  // - newPt and new segment [Barycenter-FatherNewSeg]
  if ((pNew-myVectSegments[nearIndex].myCoordinate).norm() < minDistance ||
      (pNew-myVectSegments[myVectParent[nearIndex]].myCoordinate).norm() < minDistance) {
    //DGtal::trace.warning() << "initial too close to new!!!!!!!!" << std::endl;
    return true;
  }
  // Check barycenter is not too close considered segment
  if ((pCenter-myVectSegments[nearIndex].myCoordinate).norm() < minDistance ||
      (pCenter-myVectSegments[myVectParent[nearIndex]].myCoordinate).norm() < minDistance ||
      (pCenter - pNew).norm() < minDistance){
    //DGtal::trace.warning() << "new barycenter too close!!!!!!!!" << std::endl;
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
  
  if (getProjDistanceDico(nearIndex, pNew, pCenter) < minDistance||
      getProjDistanceDico(myVectParent[nearIndex], pNew, pCenter) < minDistance||
      getProjDistanceDico(myVectChildren[nearIndex].first, pNew, pCenter) < minDistance||
      getProjDistanceDico(myVectChildren[nearIndex].second, pNew, pCenter) < minDistance){
    //DGtal::trace.warning() << "detection too close existing" << std::endl;
    return true;
  }
  
  return false;
}

bool
CoronaryArteryTree::isIntersecting(unsigned int index1, unsigned int index2, double epsilon) const
{
  Point3D p10 = myVectSegments[index1].myCoordinate;
  Point3D p11 = myVectSegments[myVectParent[index1]].myCoordinate;
  Point3D p20 = myVectSegments[index2].myCoordinate;
  Point3D p21 = myVectSegments[myVectParent[index2]].myCoordinate;
  Point3D c1 = (p10+p11)/2.0;
  Point3D c2 = (p20+p21)/2.0;
  double d = (c1-c2).norm();
  double l1 = (p10-p11).norm();
  double l2 = (p20-p21).norm();
  double tmp = d-l1/2.0-l2/2.0;
  double r1 = myVectSegments[index1].myRadius;
  double r2 = myVectSegments[index2].myRadius;
  if(tmp>r1+r2)
    return false;
  double d1 = getProjDistance(index1, p20);
  double d2 = getProjDistance(index1, p21);
  double d3 = std::min(d1,d2);
  Point3D p, q, r;
  if(d1<d2) {
    q = p21;
    p = p20;
  }
  else {
    p = p21;
    q = p20;
  }
  r = (1.0-epsilon)*p + epsilon*q;
  double d4 = getProjDistance(index1, r);
  if ((d4 > d3) && (d3 > r1+r2))
    return false;
  //Use dico search for testing intersection
  double d5 = getProjDistanceDico(index1, p20, p21, epsilon);
  return (d5 < r1+r2);
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
  myBoard.setPenColor(DGtal::Color::Blue);
  myBoard.setLineWidth(myVectSegments[0].myRadius*scaleBoard*thickness);
  myBoard.drawCircle(myTreeCenter[0], myTreeCenter[1], my_rPerf, 1);
  
  
  // draw root: (first segment is special reduced to one point and no parent).
  Point3D p0 = myVectSegments[0].myCoordinate;
  myBoard.setPenColor(DGtal::Color(10, 100, 0, 180));
  myBoard.setLineWidth(1.0);
  myBoard.fillCircle(p0[0], p0[1], myVectSegments[0].myRadius*thickness, 1);
  
  Point3D p1 = myVectSegments[1].myCoordinate;
  
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
    Point3D distal = s.myCoordinate;
    Point3D proxital = myVectSegments[myVectParent[s.myIndex]].myCoordinate;
    myBoard.setLineWidth(myVectSegments[s.myIndex].myRadius*scaleBoard*thickness);
    // myBoard.setPenColor(cmap_grad(i));
    myBoard.setPenColor(DGtal::Color::Gray);
    myBoard.drawLine(distal[0], distal[1], proxital[0], proxital[1],2);
    myBoard.setLineWidth(1.0);
    myBoard.setPenColor(DGtal::Color(180, 0, 0, 180));
    myBoard.fillCircle(distal[0], distal[1], myVectSegments[s.myIndex].myRadius*thickness, 1 );
    i++;
  }
  /*
  //Test: show constraint seeds
  for(int i=0; i<6; i++) {
    myBoard.setPenColor(DGtal::Color(0, 0, 180, 180));
    auto s = myVectSegments.at(i);
    Point3D proxital = myVectSegments[myVectParent[s.myIndex]].myCoordinate;
    myBoard.fillCircle(proxital[0], proxital[1], 2, 1);
  }
  */
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
CoronaryArteryTree::getNearestSegment(const Point3D &pt)
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
CoronaryArteryTree::getPathToRoot(const Segment<Point3D> &s)
{
  std::vector<unsigned int> res;
  Segment<Point3D> sC = s;
  while (sC.myIndex != 0){
    res.push_back(sC.myIndex);
    sC = myVectSegments[ myVectParent[sC.myIndex]];
  }
  return res;
}

CoronaryArteryTree::Point3D
CoronaryArteryTree::generateNewLocation(unsigned int nbTrials){
  Point3D res;
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


std::pair<CoronaryArteryTree::Point3D, bool>
CoronaryArteryTree::generateALocation(double myDThresold) {
  Point3D res = generateRandomPtOnDisk(myTreeCenter, my_rPerf);
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
  return  std::pair<Point3D, bool> {res, isComp};
}

bool
CoronaryArteryTree::isToCloseFromNearest(const Point3D &p, double minDist) const{
  double d = getProjDistance(getN_NearestSegments(p,1)[0],p);
  return d < minDist;
}
double
CoronaryArteryTree::getProjDistance(const Point3D &p0, const Point3D &p1, const Point3D &p) const{
  double result = 0.0;
  Point3D pProj(0,0);
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
CoronaryArteryTree::getProjDistance(unsigned int index, const Point3D &p ) const {
  Point3D p0 = myVectSegments[index].myCoordinate;
  Point3D p1 = myVectSegments[myVectParent[index]].myCoordinate;
  return  getProjDistance(p0, p1, p);
}

double
CoronaryArteryTree::getProjDistanceDico(unsigned int index1, const Point3D &p1, const Point3D &p2, const double& epsilon) const {
  Point3D p10 = myVectSegments[index1].myCoordinate;
  Point3D p11 = myVectSegments[myVectParent[index1]].myCoordinate;
  Point3D p20 = p1;
  Point3D p23 = p2;
  Point3D p21 = (2.0*p20+p23)/3.0;
  Point3D p22 = (p20+2.0*p23)/3.0;
  double d20, d21, d22, d23, delta_d1, delta_d2, delta_d3;
  double d = (p20-p22).norm();
  while (d>epsilon) {
    d20 = getProjDistance(p10, p11, p20);
    d21 = getProjDistance(p10, p11, p21);
    d22 = getProjDistance(p10, p11, p22);
    d23 = getProjDistance(p10, p11, p23);
    delta_d1 = d21 - d20;
    delta_d2 = d22 - d21;
    delta_d3 = d23 - d22;
    if(delta_d2>=0 && delta_d3>=0)
      p23 = p22;
    else
      p20 = p21;
      
    p21 = (2.0*p20+p23)/3.0;
    p22 = (p20+2.0*p23)/3.0;
    d = (p20-p22).norm();
  }
  return  std::min(d20, std::min(d21, std::min (d22, d23)));
}

std::vector<unsigned int>
CoronaryArteryTree::getN_NearestSegments(const Point3D &p, unsigned int n) const {
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
CoronaryArteryTree::hasNearestIntersections(const Point3D &p0,
                                            const Point3D &p1, unsigned int n) const {
  Point3D b = (p1+p0)/2;
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

bool CoronaryArteryTree::hasIntersections(const Point3D &p) const {
  for(auto s : myVectSegments) {
    double d = getProjDistance(s.myIndex, p);
    if (d<2*s.myRadius)
      return true;
  }
  return false;
}

bool CoronaryArteryTree::hasIntersections(unsigned int indexSeg, double epsilon) const {
  if(indexSeg==0 || indexSeg==1)
    return false;
  for(auto s : myVectSegments) {
    if (s.myIndex != 0) {
      //std::cout<<"Current: indexSeg="<<indexSeg<<" vs s.myIndex="<<s.myIndex<<std::endl;
      //std::cout<<"brother: "<<myVectSegments[myVectParent[indexSeg]].myIndex<<" vs. "<<myVectSegments[myVectParent[s.myIndex]].myIndex<<std::endl;
      //std::cout<<"children: "<<myVectSegments[myVectParent[s.myIndex]].myIndex<<" vs. "<<indexSeg<<std::endl;
      if(indexSeg != s.myIndex && // current segment
         myVectSegments[myVectParent[indexSeg]].myIndex != myVectSegments[myVectParent[s.myIndex]].myIndex && //brother = same parents
         myVectSegments[myVectParent[s.myIndex]].myIndex != indexSeg && // children
         s.myIndex != myVectSegments[myVectParent[indexSeg]].myIndex ) { //parent
        if(isIntersecting(indexSeg, s.myIndex, epsilon))
          return true;
      }
    }
  }
  return false;
}


bool
CoronaryArteryTree::hasNearestIntersections(unsigned int indexPFather,
                                            unsigned int indexPChild,
                                            const Point3D &pAdded,
                                            const Point3D &pBifurcation, unsigned int n) const{
  std::vector<unsigned int> near = getN_NearestSegments(pAdded, n);
  Point3D p0  = myVectSegments[indexPFather].myCoordinate;
  Point3D p1  = myVectSegments[indexPChild].myCoordinate;
  
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
CoronaryArteryTree::kamyiaOptimization(const Point3D& pCurrent,
                                       const Point3D& pParent,
                                       double rCurrent,
                                       const Segment<Point3D>& sL,
                                       const Segment<Point3D>& sR,
                                       unsigned int nbIter,
                                       Point3D& pOpt,
                                       double& r0, double& r1, double& r2){
  
  // Prepare
  Point3D pL = sL.myCoordinate;
  Point3D pR = sR.myCoordinate;
  
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
  //Point3D pb ((f0*pParent[0]+f1*pL[0]+f2*pR[0])/(2.0*f0), (f0*pParent[1]+f1*pL[1]+f2*pR[1])/(2.0*f0));
  //Point3D pb ((pParent[0]+pL[0])/(2.0), (pParent[1]+pL[1])/(2.0)); //center of the old segment
  Point3D pb (pCurrent[0], pCurrent[1], pCurrent[2]);
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
    hasSolution = kamiyaOpt(my_gamma, deltaP1, deltaP2, f0, f1, f2, l0, l1, l2,rr1, rr2);
    // Equation 27
    R0 = pow(f0*(pow(rr1, my_gamma)/f1 + pow(rr2, my_gamma)/f2), 1.0/my_gamma);
    R1 = rr1;
    R2 = rr2;
    // Equation (26) page 13
    pb[0] = (pParent[0]*R0/l0 + pL[0]*R1/l1 + pR[0]*R2/l2)/(R0/l0+R1/l1+R2/l2);
    pb[1] = (pParent[1]*R0/l0 + pL[1]*R1/l1 + pR[1]*R2/l2)/(R0/l0+R1/l1+R2/l2);
    pb[2] = (pParent[2]*R0/l0 + pL[2]*R1/l1 + pR[2]*R2/l2)/(R0/l0+R1/l1+R2/l2);
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
  return pow(((4.0*M_PI*my_rPerf*my_rPerf*my_rPerf)/(3.0*myKTerm)),1.0/3.0); //2D: sqrt((M_PI*my_rPerf*my_rPerf)/myKTerm);
  //Scale factor must be taken into account here
  //return sqrt((M_PI*(myKTerm + 1)*myRsupp*myRsupp)/myKTerm) / myLengthFactor;
}

CoronaryArteryTree::Point3D
CoronaryArteryTree::findBarycenter(const Point3D &p, unsigned int index)
{
  Point3D first_point=p;
  Point3D second_point=myVectSegments[index].myCoordinate;
  Point3D third_point=myVectSegments[myVectParent[index]].myCoordinate;
  Point3D barycenter((first_point[0]+second_point[0]+third_point[0])/3.0,(first_point[1]+second_point[1]+third_point[1])/3.0, (first_point[2]+second_point[2]+third_point[2])/3.0);
  return barycenter;
}

bool operator==(CoronaryArteryTree::Segment<CoronaryArteryTree::Point3D>  S1, CoronaryArteryTree::Segment<CoronaryArteryTree::Point3D> S2)
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

bool operator!=(CoronaryArteryTree::Segment<CoronaryArteryTree::Point3D>  S1, CoronaryArteryTree::Segment<CoronaryArteryTree::Point3D> S2)
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
