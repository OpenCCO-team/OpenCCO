
#include "CoronaryArteryTree.h"
#include "DGtal/io/boards/Board2D.h"
#include "geomhelpers.h"
#include "DGtal/io/readers/PPMReader.h"
#include "DGtal/images/ArrayImageIterator.h"
#include "DGtal/io/readers/PGMReader.h"
#include <math.h>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

using namespace std;
bool
CoronaryArteryTree::addSegmentFromPoint(const Point2D &p)
{
  unsigned int nearIndex = getNearestSegment(p);
  return addSegmentFromPoint(p, nearIndex);
}

bool
CoronaryArteryTree::addFirstSegment(const Point2D &p)
{
  Segment<Point2D> s;
  s.myCoordinate = p;
  s.radius = 0.1;
  s.index = myVectSegments.size();
  s.Generation=0;
  myVectSegments.push_back(s);
  myVectParent.push_back(0);
  myVectDaughters.push_back(SegmentDaugther(0,0));
    return true;
}

bool
CoronaryArteryTree::addSegmentFromPoint(const Point2D &p, unsigned int nearIndex)
{

  // To add a new segment we need to add two new points (p and middle point) with two segments.
  // (a) First point added : point associated to the nearest segment. (basic solution the center of the nearest)
  // (b) Second point added: the new point with its new segment.
  
  // to process (a): s middle
  Segment<Point2D> sMiddle;
  sMiddle.myCoordinate = FindBarycenter(p, nearIndex);
  sMiddle.radius = 0.1;
//  sMiddle.index = myVectSegments.size(); modifié
//  myVectSegments.push_back(sMiddle); modifié
//  myVectParent.push_back(myVectParent[nearIndex]); modifié
    
    Segment<Point2D> sNew2;
    sMiddle.index = nearIndex;
    sNew2 = myVectSegments[nearIndex];
    sNew2.index = myVectSegments.size();
    myVectSegments.push_back(sNew2);
    myVectParent.push_back(nearIndex);
    myVectSegments[nearIndex] = sMiddle;
    
  
    myVectDaughters.push_back(SegmentDaugther(sNew2.index,  sNew2.index)); // modifié

  // to process (b): new point to s middle
  Segment<Point2D> sNew;
  sNew.myCoordinate = p;
  sNew.radius = 0.1;
  sNew.index = myVectSegments.size();
  myVectSegments.push_back(sNew);
  myVectParent.push_back(sMiddle.index);// Modifié
    
  // update parent with new (a).
 //  myVectParent[nearIndex] = sMiddle.index; modifié
    //We define the Daughter of a Terminal segment as himself
    myVectDaughters.push_back(SegmentDaugther(sNew.index, sNew.index));
    myVectDaughters[nearIndex] = SegmentDaugther( sNew2.index, sNew.index);
  // NearIndex est un fils, sNew est un fils, et sMiddle est le pere des deux.
    updateGeneration();
    updateRadius2(sNew.index);
  
  return true;
}



void
CoronaryArteryTree::exportDisplay(const std::string &fileName)
{
  // drawing base circle
  DGtal::Board2D board;
    std::cout <<"My rsupp is " << myRsupp<<std::endl;
  board.setPenColor(DGtal::Color::Blue);
  board.drawCircle(myTreeCenter[0], myTreeCenter[1], myRsupp);
  board.setPenColor(DGtal::Color::Red);
  

  // draw root: (first segment is special reduced to one point and no parent).
  Point2D p0 = myVectSegments[0].myCoordinate;
  board.setPenColor(DGtal::Color::Black);
//  board.setFillColor(DGtal::Color::Red);
  board.fillCircle(p0[0], p0[1], 0.5);
  
  Point2D p1 = myVectSegments[1].myCoordinate;
  board.setLineWidth(myVectSegments[0].radius*3.0);
//  board.setPenColor(DGtal::Color::Red);
  board.drawLine(p0[0], p0[1], p1[0], p1[1]);
  
  for (auto s : myVectSegments)
  {
    // test if the segment is the root or its parent we do not display (already done).
    if (s.index == 0 || s.index == 1)
    {
      continue;
    }
    // distal node:
    Point2D distal = s.myCoordinate;
    Point2D proxital = myVectSegments[myVectParent[s.index]].myCoordinate;
    board.setPenColor(DGtal::Color::Black);
    //board.setFillColor(DGtal::Color::Red);
    //board.fillCircle(distal[0], distal[1], 0.5);
      board.setLineWidth(sqrt(s.radius*50));
    board.setPenColor(DGtal::Color::Gray);
    board.drawLine(distal[0], distal[1], proxital[0], proxital[1]);
//      std::cout << " myRsupp Value = "<< myRsupp<<std::endl;
//      board.drawCircle(distal[0], distal[1],myDThresold);

    
  }
  
  std::string ext = fileName.substr(fileName.find_last_of(".") + 1);
  if (ext == "svg")
  {
    board.saveSVG(fileName.c_str());
  }
  else if (ext == "eps") {
    board.saveEPS(fileName.c_str());
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
  if (s.index != 0)
  {
    unsigned int p = myVectParent.at(s.index);
    res += myVectSegments.at(p).myCoordinate;
    res[0] /= 2.0;
    res[1] /= 2.0;
  }
  return res;
}


unsigned int
CoronaryArteryTree::getNearestSegment(const Point2D &pt)
{
  unsigned int sNear=0;
  for (const auto &s: myVectSegments)
  {
    auto c = getSegmentCenter(s);
    auto distMin = (getSegmentCenter(myVectSegments[sNear])-pt).norm();
    if ((getSegmentCenter(s)-pt).norm() < distMin)
    {
      sNear = s.index;
    }
  }
  return sNear;
}



unsigned int
CoronaryArteryTree::getParentSegment(const Segment<Point2D> &s){
  return myVectParent[s.index];
}

unsigned int
CoronaryArteryTree::getDaughterLeft(const Segment<Point2D> &s){
  return myVectDaughters[s.index].first;
}

unsigned int
CoronaryArteryTree::getDaughterRigth(const Segment<Point2D> &s){
  return myVectDaughters[s.index].second;
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

    Point2D first_point=p;
    Point2D second_point=myVectSegments[nearIndex].myCoordinate;
    Point2D third_point=myVectSegments[myVectParent[nearIndex]].myCoordinate;
    Point2D barycenter=FindBarycenter(p,nearIndex);
    
    Segment<Point2D> sMiddle;
    sMiddle.myCoordinate = barycenter;
    sMiddle.radius = sqrt(qTerm/(pi*GetLength(nearIndex)));
    sMiddle.index = myVectSegments.size();
    myVectSegments.push_back(sMiddle);
    myVectParent.push_back(myVectParent[nearIndex]);
    myVectDaughters.push_back(SegmentDaugther(nearIndex,  myVectSegments.size()));

    // to process (b): new point to s middle
    Segment<Point2D> sNew;
    sNew.myCoordinate = p;
    sNew.radius = sqrt(qTerm/(pi*GetLength(nearIndex)));
    sNew.index = myVectSegments.size();
    myVectSegments.push_back(sNew);
    myVectParent.push_back(sMiddle.index);
    // update parent with new (a).
    myVectParent[nearIndex] = sMiddle.index;
    
    // NearIndex est un fils, sNew est un fils, et sMiddle est le pere des deux.
    myVectDaughters[sMiddle.index] = SegmentDaugther(nearIndex, sNew.index);
    
    updateGeneration();
    
    updateRadius();

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

bool
CoronaryArteryTree::updateRsupp()
{
    myRsupp=sqrt(myKTerm*(aPerf/Nterm)/pi);
//    cout << "My new Rsupp is = " << myRsupp << endl;
    return true;
}

bool
CoronaryArteryTree::updateTreshold()
{
    myDThresold=sqrt(pi*myRsupp*myRsupp/myKTerm);
    return true;
}
bool
CoronaryArteryTree::updateRadius()
{
    for (auto s : myVectSegments)
    {
        
        myVectSegments[s.index].radius=sqrt((1+MaxGene-s.Generation)*qTerm/(pi*GetLength(s.index)));
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


        myVectSegments[getDaughterLeft(parent)].radius = sqrt((1+MaxGene- myVectSegments[getDaughterLeft(parent)].Generation)*qTerm/(pi*GetLength(getDaughterLeft(parent))));
    
    myVectSegments[getDaughterRigth(parent)].radius = sqrt((1+MaxGene-myVectSegments[getDaughterRigth(parent)].Generation)*qTerm/(pi*GetLength(getDaughterRigth(parent))));
    // Now, we need to change the parents according to the Bifurcation rule;
    
    while(parent.Generation != 0)
    {
    myVectSegments[parent.index].radius = pow(pow(myVectSegments[getDaughterRigth(parent)].radius,gamma) + pow(myVectSegments[getDaughterLeft(parent)].radius,gamma),1/gamma);
        
        parent=myVectSegments[myVectParent[parent.index]];

    }
    myVectSegments[parent.index].radius =sqrt((myKTerm)*qTerm/(pi*GetLength(parent.index)));
    
    return TRUE;

}




bool
CoronaryArteryTree::updateScale( double rescale_factor)
{
    updateRsupp();
    
    updateTreshold();
    for (auto s : myVectSegments)
    {
        s.myCoordinate*=(myRsupp/rescale_factor);

    }
    
    return true;
}

bool
CoronaryArteryTree::updateGeneration()

{
    int Gene=1;
    Segment<Point2D> parent;
    for (auto s : myVectSegments)
    {
        if(s.index != 0){
            
        
        parent = myVectSegments[myVectParent[s.index]];
        int Gene=1;
        
        while(GetLength(parent.index)!=0 )
        {
            parent = myVectSegments[myVectParent[parent.index]];
            Gene+=1;
        }
        myVectSegments[s.index].Generation = Gene;
        if(Gene>MaxGene){MaxGene=Gene;}
//        s.Generation = Gene;
        
    }
    }
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


CoronaryArteryTree::Point2D
CoronaryArteryTree::FindOptimalPosition(unsigned int Index,const Point2D &p)
{
    // We have to find the best postion, without drawing any new segment.
    
    //First, we have 3 points: p, the end of the Index segment, and the end of his father
    
    double totalVolume;
    Point2D first_point=p; // This is the point to add
    Point2D second_point=myVectSegments[Index].myCoordinate; // This is the future parent of the new point
    Point2D third_point=myVectSegments[myVectParent[Index]].myCoordinate;// This is the parent of the segment on which we want to connect the point

    Point2D OptimalPosition= getSegmentCenter(Index); // We initialiaze the gradient descent with this point, it's the segment center.

    double H1,H2,H3;
    H1= sqrt(pow((first_point[0]-OptimalPosition[0]),2) + pow(first_point[1]-OptimalPosition[1],2));
    H2= sqrt(pow((second_point[0]-OptimalPosition[0]),2) + pow(second_point[1]-OptimalPosition[1],2));
    H3= sqrt(pow((third_point[0]-OptimalPosition[0]),2) + pow(third_point[1]-OptimalPosition[1],2));
    
//    totalVolume = pi*myVectSegments[Index].radius*myVectSegments[Index].radius*H +
//    pi*myVectSegments[myVectParent[nearIndex]].radius*myVectSegments[myVectParent[nearIndex]].radius*H
//    +pi*myVectSegments[Index].radius*myVectSegments[Index].radius*H
    totalVolume = pi*H1 + pi*H2+pi*H3;
    //
    double xnew;
    double ynew;
    xnew=OptimalPosition[0];
    ynew=OptimalPosition[1];
    double VolumeXDerivative=pi*(2*(xnew - first_point[0])* (1/(2*(sqrt(pow((first_point[0]-xnew),2) +pow(first_point[1]-ynew,2))))))
    + pi*(2*(xnew - second_point[0])* (1/(2*(sqrt(pow((second_point[0]-xnew),2) + pow(second_point[1]-ynew,2))))))
    + pi*(2*(xnew - third_point[0])* (1/(2*(sqrt(pow((third_point[0]-xnew),2) + pow(third_point[1]-ynew,2))))));
    
    double VolumeYDerivative=pi*(2*(ynew - first_point[1])* (1/(2*(sqrt(pow((first_point[0]-xnew),2) +pow(first_point[1]-ynew,2))))))
    + pi*(2*(ynew - second_point[1])* (1/(2*(sqrt(pow((second_point[0]-xnew),2) + pow(second_point[1]-ynew,2))))))
    + pi*(2*(ynew - third_point[1])* (1/(2*(sqrt(pow((third_point[0]-xnew),2) + pow(third_point[1]-ynew,2))))));

    
    // We have the value of the derivative in X and Y
    double alpha=0.0001;
    int nb_itarations=0;
    double grad = 1000;
    double grad2 = 1000;

    //&& nb_itarations<1000
//    std::cout << "bool value is " <<(abs(grad-abs(pow(pow(VolumeYDerivative,2)+pow(VolumeXDerivative,2),0.5)) )<0.001) <<std::endl;
    while(abs(pow(pow(VolumeYDerivative,2)+pow(VolumeXDerivative,2),0.5)) > 0.5 && abs(grad -abs(pow(pow(VolumeYDerivative,2)+pow(VolumeXDerivative,2),0.5)))>0.0001 ){
    nb_itarations++;
    double dxH1,dxH2,dxH3,dyH1,dyH2,dyH3;
    H1= sqrt(pow((first_point[0]-xnew),2) + pow(first_point[1]-ynew,2));
    H2= sqrt(pow((second_point[0]-xnew),2) + pow(second_point[1]-ynew,2));
    H3= sqrt(pow((third_point[0]-xnew),2) + pow(third_point[1]-ynew,2));
        
    
        
    xnew= xnew -  alpha*VolumeXDerivative;
    ynew= ynew - alpha*VolumeYDerivative;
    

    totalVolume = pi*H1 + pi*H2+pi*H3;
    grad2=abs(pow(pow(VolumeYDerivative,2)+pow(VolumeXDerivative,2),0.5));
        
}
    Point2D newpoint(xnew,ynew);
    std::cout << " il y a eu " << nb_itarations << " iterations "<<std::endl;
    return newpoint;
}

double
CoronaryArteryTree::GetTotalVolume(const Point2D &p1,const Point2D &p2,const Point2D &p3,const Point2D &pOpti)
{
    double H1,H2,H3;
    
    H1 = sqrt(pow(p1[0]-pOpti[0],2)+pow(p1[1]-pOpti[1],2));
    H2 = sqrt(pow(p2[0]-pOpti[0],2)+pow(p2[1]-pOpti[1],2));
    H3 = sqrt(pow(p3[0]-pOpti[0],2)+pow(p3[1]-pOpti[1],2));
    
    return pi*H1+pi*H2+pi*H3;
    
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
  sMiddle.radius = 0.1;
  sMiddle.index = myVectSegments.size();
  myVectSegments.push_back(sMiddle);
  myVectParent.push_back(myVectParent[nearIndex]);
  myVectDaughters.push_back(SegmentDaugther(nearIndex,  myVectSegments.size()));

  // to process (b): new point to Optimal position
  Segment<Point2D> sNew;
  sNew.myCoordinate = NewPoint;
  sNew.radius = 0.1;
  sNew.index = myVectSegments.size();
  myVectSegments.push_back(sNew);
  myVectParent.push_back(sMiddle.index);
  // update parent with new (a).
  myVectParent[nearIndex] = sMiddle.index;


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
CoronaryArteryTree::FindOptimalSegment(const Point2D &p)
{
    double Volume_min=1000000;
    int Index_Opti;
    DGtal::Z2i::RealPoint FinalPosition;
    
    for (auto s : myVectSegments)
    {
        DGtal::Z2i::RealPoint OptimalPosition;
        OptimalPosition = FindBarycenter(p,s.index);
        
        if(compDistCriteria(p, s.index) > myDThresold){
            
        
        if(Volume_min>GetTotalVolume( s.myCoordinate, myVectSegments[myVectParent[s.index]].myCoordinate, p,OptimalPosition )&& s.index >0 )
        {
            
            Volume_min =GetTotalVolume(s.myCoordinate, myVectSegments[myVectParent[s.index]].myCoordinate, p,OptimalPosition );
            Index_Opti= s.index;
            
        
        }
        }
        else
        {
            return -1;
        }

    }
//    cout << " we return " << Index_Opti<<endl;
    return Index_Opti;
}

int
CoronaryArteryTree::AddFirstSegmentonImage()
{
    DGtal::Z2i::RealPoint p(320,200);

    Segment<Point2D> s;
    s.myCoordinate = p;
    s.radius = 0.1;
    s.index = myVectSegments.size();
    s.Generation=1;
    myVectSegments.push_back(s);
    myVectParent.push_back(0);
    myVectDaughters.push_back(SegmentDaugther(0,0));
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

//
//double
//CoronaryArteryTree::FindEdge(int xDim, int yDim)
//{
//
//    typedef DGtal::ImageSelector < DGtal::Z2i::Domain, unsigned int>::Type Image;
//    Image image = DGtal::PGMReader<Image>::importPGM( "./Shape.pgm" );
//
//    for(auto const &point: image.domain())
//    {
//        if(point == 0)
//        {
//
//            int ximage;
//            int yimage;
//            ximage = -(x-xmax)*(image.domain().upperBound()[0]/2)/xmax;
//            yimage = -(y-ymax)*(image.domain().upperBound()[1]/2)/ymax;
//            return
//        }
//    }
//    return myRsupp*sin(atan(xDim/yDim));
//
//}
//

// The following function converts a pixelfrom the image into a point in the supporting circle
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



// The following function converts a random generated point into a couple (x,y) in the picture
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
   DGtal::Z2i::RealPoint Barycenter = FindBarycenter(newPoint, S1.index);
   
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
           
       
           if(hasIntersection(newPoint ,Barycenter,s.myCoordinate,myVectSegments[myVectParent[s.index]].myCoordinate) || hasIntersection(s.myCoordinate,myVectSegments[myVectParent[s.index]].myCoordinate,Barycenter,S1.myCoordinate) || hasIntersection(s.myCoordinate,myVectSegments[myVectParent[s.index]].myCoordinate,Barycenter,myVectSegments[myVectParent[S1.index]].myCoordinate) )
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
