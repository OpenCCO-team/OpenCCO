//
//  testResize.cpp
//  
//
//  Created by Dylan on 24/06/2020.
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "CoronaryArteryTree.h"

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/boards/Board2D.h"
#define PI 3.14159265

#include "geomhelpers.h"
using namespace DGtal;
using namespace DGtal::Z2i;
using namespace std;

int main(int argc, char *const *argv)
{
    int nb_segment=5;
  bool ok = true;
  DGtal::trace.beginBlock("Testing class CoronaryArteryTree: (Resize)");
  CoronaryArteryTree c (DGtal::Z2i::RealPoint(0, 25), 2000, 1);
  DGtal::trace.info() << c;
  c.Nterm=nb_segment;
  
//    cout << " myRsupp value = " << c.myRsupp <<endl;
//    cout << " myKTerm value = " << c.myKTerm <<endl;
//    cout << " myDTresholsd = " << c.myDThresold <<endl;
//    c.myKTerm++;
//    c.updateRsupp();
//    c.updateTreshold();
//    cout << " after update, myRsupp value = " << c.myRsupp <<endl;
//    cout << " after update, myKTerm value = " << c.myKTerm <<endl;
//    cout << " myDTresholsd = " << c.myDThresold <<endl;
//
//    c.myKTerm++;
//    c.updateRsupp();
//    c.updateTreshold();
//    cout << " after update, myRsupp value = " << c.myRsupp <<endl;
//    cout << " after update, myKTerm value = " << c.myKTerm <<endl;
//    cout << " myDTresholsd = " << c.myDThresold <<endl;
//
//    board.setLineWidth(3.0);
//    board.setPenColor(DGtal::Color::Red);
//
//      board.drawCircle(0, 0, c.myRsupp);
      
      c.exportDisplay("testStep2Resize.svg");
    

    srand (time(NULL));
    CoronaryArteryTree cRand (DGtal::Z2i::RealPoint(0, 25), 2000, 1);
    cRand.Nterm=200;
    cRand.addFirstSegment(DGtal::Z2i::RealPoint(0, 0));
    for (unsigned int i = 0; i < 200; i++){
      double r=( rand()/double(RAND_MAX))*cRand.myRsupp;
      double theta =(rand()/double(RAND_MAX))*2*PI;
          
      double x= r*cos(theta);
      double y= r*sin(theta);
      if (isInsideCircle(DGtal::Z2i::RealPoint(0, 0), DGtal::Z2i::RealPoint(x, y),  cRand.myRsupp)){
        double oldRsupp=cRand.myRsupp;
        cRand.addSegmentFromPoint(DGtal::Z2i::RealPoint(x, y));
        cRand.myKTerm++;
//        cRand.updateRsupp();
//        cRand.updateTreshold();
          
        cRand.updateScale(oldRsupp);
          
      }
    }
    
    
     //Now, let's try to add a new point with barycenter method
    bool ok5=true;
    double x,y,r,theta;
    while(ok5){
    DGtal::trace.beginBlock("Testing class CoronaryArteryTree: test gradient descent ");
        r=( rand()/double(RAND_MAX))*cRand.myRsupp;
        theta =(rand()/double(RAND_MAX))*2*PI;
           
        x= r*cos(theta);
        y= r*sin(theta);
        if (isInsideCircle(DGtal::Z2i::RealPoint(0, 0), DGtal::Z2i::RealPoint(x, y),  cRand.myRsupp))
        {ok5=false;}
        
            }
    DGtal::Z2i::RealPoint NewPoint(x,y);
    cout << "le nouveau point est : x = " << x << " and y = " << y <<endl;
    std::vector<unsigned int >  Elligible_Index;
    cout << cRand.myDThresold <<endl;
    //Parcourons l'arbre pour trouver a quel branche s'attacher
    for(int i=0;i<cRand.myVectSegments.size();i++){
    if(cRand.compDistCriteria(DGtal::Z2i::RealPoint(x, y), i)>cRand.myDThresold){
            Elligible_Index.push_back(i);
        }
    }
    cout << "elligible size = " << Elligible_Index.size()<<endl;
    
    if(Elligible_Index.size()>0){
        double Volume_min=100000;
        int Index_Opti;
        bool test = false;
        DGtal::Z2i::RealPoint FinalPosition;
    for(int i=0;i<Elligible_Index.size();i++)
    {
        DGtal::Z2i::RealPoint OptimalPosition;
        OptimalPosition = cRand.FindBarycenter(NewPoint,i );
        if(Volume_min>cRand.GetTotalVolume( cRand.myVectSegments[Elligible_Index[i]].myCoordinate, cRand.myVectSegments[cRand.myVectParent[Elligible_Index[i]]].myCoordinate, NewPoint,OptimalPosition ))
                {
                    cout << "UPDATE ! New volume value = " <<cRand.GetTotalVolume( cRand.myVectSegments[Elligible_Index[i]].myCoordinate, cRand.myVectSegments[cRand.myVectParent[Elligible_Index[i]]].myCoordinate, NewPoint,OptimalPosition )<<endl;
                  
                    Volume_min =cRand.GetTotalVolume( cRand.myVectSegments[Elligible_Index[i]].myCoordinate, cRand.myVectSegments[cRand.myVectParent[Elligible_Index[i]]].myCoordinate, NewPoint,OptimalPosition );
                    cout << " the index is " << Index_Opti << " and the point is x = "<< OptimalPosition[0] << " and y = " << OptimalPosition[1] << endl;
                    Index_Opti=i;
                    FinalPosition=OptimalPosition;

                    test=TRUE;
                }
    
    
    }
        
    }
    
    
    cout << "Kterm = " << cRand.myKTerm<<endl;
    cout << "rsupp = " << cRand.myRsupp<<endl;

    cout << "Thresold = " << cRand.myDThresold<<endl;
    cRand.exportDisplay("testResize.svg");
    
    
    return 0;
}


  
    
  


