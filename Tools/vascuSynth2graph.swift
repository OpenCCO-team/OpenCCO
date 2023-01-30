//
//  main.swift
//  vascuSynth2graph
//
//  Created by Bertrand Kerautret on 10/06/2022.
//

import Foundation

struct Vertex {
  var id: String? = nil
  var x: Double? = nil
  var y: Double? = nil
  var z: Double? = nil
  var radius: Double? = nil
}

struct Edge {
  var idA: String?
  var idb: String?
  var radius: Double?
}


class FileVascuParse : NSObject, XMLParserDelegate {
  //var description: String
  var myURL: URL
  
  var myVertices = [Vertex]()
  var myEdges = [Edge]()
  var myCurrentVertex: Vertex?
  var myCurrentEdge: Edge?
  var inNode = false
  var readRadius = false
  var nodeId: String?
  var foundCharacters = ""
  init(anURL: URL) {
    myURL = anURL
  }
  func apply(){
    let parser = XMLParser(contentsOf: myURL)!
    parser.delegate = self
    let success = parser.parse()
    if success {
      print("Success...")
    }
  }
  
  func parser(_ parser: XMLParser, foundCharacters string: String) {
     
    self.foundCharacters += string
  }
  
  func parser(_ parser: XMLParser, didStartElement elementName: String, namespaceURI: String?, qualifiedName qName: String?, attributes attributeDict: [String : String] = [:]) {
    if myCurrentEdge != nil {
      if elementName == "attr" && attributeDict["name"] == " radius"{
        readRadius = true
      }
    }
    if elementName == "node" {
      myCurrentVertex = Vertex()
      inNode = true
      nodeId = attributeDict["id"]! as String
      myCurrentVertex?.id = nodeId
    }
    if elementName == "edge" {
      myCurrentEdge = Edge()
      myCurrentEdge?.idA = attributeDict["from"]!
      myCurrentEdge?.idb = attributeDict["to"]!
    }
      

  }
  
  func parser(_ parser: XMLParser, didEndElement elementName: String, namespaceURI: String?, qualifiedName qName: String?) {
   
    if elementName == "node" {
      if myCurrentVertex != nil && myCurrentVertex!.x != nil &&
          myCurrentVertex!.y != nil && myCurrentVertex!.z != nil &&
          myCurrentVertex?.id != nil
      {
        inNode = false
          if nodeId == "n0"{
              print("root index: \(myEdges.count)")
          }
        myVertices.append(myCurrentVertex!)
        print("item addedv \(myCurrentVertex!.id!)")
        myCurrentVertex = nil
      }
       
    }
    if elementName == "edge" {
      myEdges.append(myCurrentEdge!)
      myCurrentEdge = nil
    }
    if elementName == "float" && readRadius {
      foundCharacters = String(foundCharacters.filter { !" \n\t\r".contains($0) })

      myCurrentEdge?.radius = Double(foundCharacters)
      foundCharacters=""
      readRadius = false
    }
    if elementName == "float" && myCurrentVertex != nil{
      foundCharacters = String(foundCharacters.filter { !" \n\t\r".contains($0) })
      if myCurrentVertex?.x == nil
      {
        myCurrentVertex?.x = Double( foundCharacters)
      }
      else if myCurrentVertex?.y == nil
      {
        myCurrentVertex?.y = Double(foundCharacters)
      }
      else if myCurrentVertex?.z == nil
      {
        myCurrentVertex?.z = Double(foundCharacters)
      }
    }
    foundCharacters=""
  }
  
}


let a = CommandLine.arguments[1]
