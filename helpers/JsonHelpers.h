
#pragma once

#if defined(JSONHELPERS_RECURSES)
#error Recursive header files inclusion detected in JSONHelpers.h
#else // defined(JSONHELPERS_RECURSES)
/** Prevents recursive inclusion of headers. */
#define JSONHELPERS_RECURSES

#if !defined JSONHELPERS_h
/** Prevents repeated inclusion of headers. */
#define JSONHELPERS_h
#include <string.h>
#include <iostream>
#include "CoronaryArteryTree.h"

using namespace std;
  
namespace JSONHelpers {

/**
 * prints a node into Json format
 */
template <int TDim>
void subPrint_node(int nodeType, int idNode, DGtal::PointVector<TDim, double> pos, double r, std::vector<unsigned int> prec, std::vector<unsigned int> next, ofstream &os){
  os<<"\t\t{\n"<<endl;
  
  os<<"\t\t\t\"id\" : "<<idNode<<","<<endl;
  os<<"\t\t\t\"centerline\" : [ ]"<<","<<endl;
  
  os<<"\t\t\t\"position\" : [ ";
  for (auto i=0; i < TDim - 1; i++){
   os<<pos[i]<<", ";
  }
  os<<pos[TDim - 1]<<"],"<<endl;
  
  os<<"\t\t\t\"prec\" : [ ";
  if(prec.size()>0) {
    for (auto i=0; i < prec.size() - 1; i++){
      os<<prec[i]<<", ";
    }
    os<<prec[prec.size() - 1]<<"],"<<endl;
  }
  else {
    os<<" ],"<<endl;
  }
  os<<"\t\t\t\"next\" : [ ";
  if(next.size()>0){
    for (auto i=0; i < next.size() - 1; i++){
      os<<next[i]<<", ";
    }
    os<<prec[next.size() - 1]<<"],"<<endl;
  }
  else {
    os<<" ],"<<endl;
  }
  
  os<<"\t\t\t\"autre_prop\" : ";
  if(nodeType == NodeTable::ROOT){
    os<<"\"root node\""<<endl;
  } else if(nodeType == NodeTable::TERM){
    os<<"\"terminal node\""<<endl;
  } else if(nodeType == NodeTable::BIF){
    os<<"\"bifurication\""<<endl;
  } else {
    os<<"\"unknown type\""<<endl;
  }
  os<<"\t\t},"<<endl;
}

/**
 * prints an edge into Json format
 */
 void subPrint_edge(int idSeg, int idSegPar, double flow, double radius, ofstream &os){
 
  if(idSeg != 0){
    os<<"\t\t{\n"<<endl;
    os<<"\t\t\t\"id\" : "<<idSeg<<","<<endl;
    os<<"\t\t\t\"from\" : "<<idSeg<<","<<endl;
    os<<"\t\t\t\"to\" : "<<idSegPar<<","<<endl;
    os<<"\t\t\t\"flow\" : "<<flow<<","<<endl;
    os<<"\t\t\t\"radius\" : "<<radius<<","<<endl;
    os<<"\t\t},"<<endl;
      
  }
}

template<int TDim>
inline
void writeTreeToJson(const CoronaryArteryTree<TDim>& tree, const char * filePath) {
  ofstream output;
  
  //writing the tree structure as GXL to the filePath specified
  output.open(filePath);
  output<<"{\n\t\"control_points\":\n\t["<<endl;
  
  //writing tree's nodes
   for(auto s : tree.myVectSegments) {
     // test if the segment is the root or its parent
     if (s.myIndex == 0) { //root node
       std::vector<unsigned int> prec = {0};
       //std::vector<unsigned int> next = {tree.myVectChildren[s.myIndex].first, tree.myVectChildren[s.myIndex].second};
       subPrint_node<TDim>(NodeTable::ROOT, s.myIndex, s.myCoordinate, s.myRadius, prec, std::vector<unsigned int>(), output);
     }
     else {
       if(std::find(begin(tree.myVectTerminals), end(tree.myVectTerminals), s.myIndex) != end(tree.myVectTerminals)) { //terminal node
         std::vector<unsigned int> prec = {tree.myVectParent[s.myIndex]};
         //std::vector<unsigned int> next = {tree.myVectChildren[s.myIndex].first, tree.myVectChildren[s.myIndex].second};
         subPrint_node<TDim>(NodeTable::TERM, s.myIndex, s.myCoordinate, s.myRadius, prec, std::vector<unsigned int>(), output);
       }
       else { // bif node
         std::vector<unsigned int> prec = {tree.myVectParent[s.myIndex]};
         std::vector<unsigned int> next = {tree.myVectChildren[s.myIndex].first, tree.myVectChildren[s.myIndex].second};
         subPrint_node<TDim>(NodeTable::BIF, s.myIndex, s.myCoordinate, s.myRadius, prec, next, output);
       }
     }
   }
  //output<<"\t]\n}, "<<endl;
  output<<"\t], "<<endl;
  
  //writing tree's edges
  //output<<"{\n\t\"edges\":\n\t["<<endl;
  output<<"\t\"edges\":\n\t["<<endl;
  for(auto s : tree.myVectSegments) {
    subPrint_edge(s.myIndex,tree.myVectParent[s.myIndex], s.myFlow, s.myRadius, output);
  }
  output<<"\t]\n}"<<endl;
  output.close();
}
};

#endif // !defined JSONHELPERS_h

#undef JSONHELPERS_RECURSES
#endif // else defined(JSONHELPERS_RECURSES)
