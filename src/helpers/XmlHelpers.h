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

#if defined(XMLHELPERS_RECURSES)
#error Recursive header files inclusion detected in Xmlhelpers.h
#else // defined(XMLHELPERS_RECURSES)
/** Prevents recursive inclusion of headers. */
#define XMLHELPERS_RECURSES

#if !defined XMLHELPERS_h
/** Prevents repeated inclusion of headers. */
#define XMLHELPERS_h
#include <string.h>
#include <iostream>
#include "CoronaryArteryTree.h"
using namespace std;

class NodeTable {
public:
    enum Possibilite {TERM, ROOT, BIF};
    /*
     //constants used for the type field of a node
     static int TERM;
     static int ROOT;
     static int BIF;
     static int FIELDS;
     
     NodeTable();
     */
};

namespace XMLHelpers {
/**
 * prints a node into XML/GXL format from the NodeTable
 */
template <typename TPoint>
void subPrint_node(int nodeType, TPoint pos,
                   int idNode, ofstream &os){
    os<<"  <node id=\"n"<<idNode<<"\">"<<endl;
    os<<"    <attr name=\"nodeType\">"<<endl;
    if(nodeType == NodeTable::ROOT){
        os<<"      <string> root node </string>"<<endl;
    } else if(nodeType == NodeTable::TERM){
        os<<"      <string> terminal node </string>"<<endl;
    } else if(nodeType == NodeTable::BIF){
        os<<"      <string> bifurication </string>"<<endl;
    } else {
        os<<"      <string> unknown type </string>"<<endl;
    }
    os<<"    </attr>"<<endl;
    
    os<<"    <attr name=\"position\">"<<endl;
    os<<"      <tup>"<<endl;
    
    for (auto i=0; i <pos.dimension; i++){
        os<<"        <float>"<<pos[i]<<"</float>"<<endl;
    }
    
    os<<"      </tup>"<<endl;
    os<<"    </attr>"<<endl;
    os<<"  </node>"<<endl;
}

/**
 * prints an edge into XML/GXL format from a node table
 */

void subPrint_edge(int idSeg, int idSegPar,
                   double flow, double radius,double resist,
                   ofstream &os){
    
    if(idSeg != 0){
        os<<"  <edge id=\"e"<<idSeg<<"\" to=\"n"<<idSeg<<"\" from=\"n"<<idSegPar<<"\">"<<endl;
        os<<"    <attr name=\"flow\">"<<endl;
        os<<"      <float>"<<flow<<"</float>"<<endl;
        os<<"    </attr>"<<endl;
        os<<"    <attr name=\"resistance\">"<<endl;
        os<<"      <float>"<<resist<<"</float>"<<endl;
        os<<"    </attr>"<<endl;
        
        os<<"    <attr name=\"radius\">"<<endl;
        os<<"      <float>"<<radius<<"</float>"<<endl;
        os<<"    </attr>"<<endl;
        
        os<<"  </edge>"<<endl;
    }
}


template<typename TTree>
inline
void writeTreeToXml(const TTree &tree, const char * filePath) {
    ofstream output;
    //writing the tree structure as GXL to the filePath specified
    output.open(filePath);
    output<<"<gxl><graph id=\""<<filePath<<"\" edgeids=\" true\" edgemode=\" directed\" hypergraph=\" false\">"<<endl;
    output<<"<info_graph>"<< endl;
    output<<"    <attr name=\"pPerf\">"<<endl;
    output<<"      <float>"<<tree.my_pPerf<<"</float>"<<endl;
    output<<"    </attr>"<<endl;
    output<<"    <attr name=\"pTerm\">"<<endl;
    output<<"      <float>"<<tree.my_pTerm<<"</float>"<<endl;
    output<<"    </attr>"<<endl;
    output<<"</info_graph>"<< endl;
    
    //writing tree's nodes
    for(auto s : tree.myVectSegments) {
        // test if the segment is the root or its parent
        if (tree.myVectParent[s.myIndex]==0) //root node
            subPrint_node(NodeTable::ROOT, s.myCoordinate, s.myIndex, output);
        else {
            if(std::find(begin(tree.myVectTerminals), end(tree.myVectTerminals), s.myIndex) != end(tree.myVectTerminals)) { //terminal node
                subPrint_node(NodeTable::TERM, s.myCoordinate, s.myIndex, output);
            }
            else { // bif node
                subPrint_node(NodeTable::BIF, s.myCoordinate, s.myIndex, output);
            }
        }
    }
    //writing tree's edges
    for(auto s : tree.myVectSegments) {
        subPrint_edge(s.myIndex,tree.myVectParent[s.myIndex], s.myFlow, s.myRadius, s.myResistance, output);
    }
    
    output<<"</graph></gxl>"<<endl;
    output.close();
}

}; // end namespace XMLHelpers

#endif // !defined XMLHELPERS_h

#undef XMLHELPERS_RECURSES
#endif // else defined(XMLHELPERS_RECURSES)
