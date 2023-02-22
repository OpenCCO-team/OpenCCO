//
//  main.swift
//  vascuSynth2graph
//
//  Created by Bertrand Kerautret on 10/06/2022.
//



import Foundation
#if canImport(FoundationXML)
 import FoundationXML
#endif

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
    var flow: Double?
}


class XMLGraphParse : NSObject, XMLParserDelegate {
    //var description: String
    var myURL: URL
    
    var myVertices = [Vertex]()
    var myEdges = [Edge]()
    var myCurrentVertex: Vertex?
    var myCurrentEdge: Edge?
    var inNode = false
    var readRadius = false
    var readFlow = false
    var nodeId: String?
    var foundCharacters = ""
    var myVerbose: Bool
    init(anURL: URL, verbose: Bool = false) {
        myURL = anURL
        myVerbose = verbose
    }
    func apply(){
        let parser = XMLParser(contentsOf: myURL)!
        parser.delegate = self
        let success = parser.parse()
        if success && myVerbose {
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
            if elementName == "attr" && attributeDict["name"] == " flow"{
                readFlow = true
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
                myCurrentVertex!.y != nil &&
                myCurrentVertex?.id != nil
            {
                inNode = false
                if nodeId == "n0"  && myVerbose {
                    print("root index: \(myEdges.count)")
                }
                myVertices.append(myCurrentVertex!)
                if myVerbose {
                    print("item addedv \(myCurrentVertex!.id!)")
                }
                myCurrentVertex = nil
            }
        }
        if readFlow && elementName == "float" {
            foundCharacters = String(foundCharacters.filter { !" \n\t\r".contains($0) })
            myCurrentEdge?.flow = Double(foundCharacters)
            foundCharacters=""
            readFlow = false
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

func usage(){
    print("Command line args: \(CommandLine.arguments[0]) input xml [outBaseName] [verbose]")
    print("Extract the xml file into basic graph representation.",
          "The representation is given in 3 files outBaseName_vertex.dat outBaseName_edges.dat outputBaseName_radius.dat")
}


if CommandLine.arguments.count < 1 {
    usage()
    exit(1)
}
var outBaseName = ""
var verbose = false
let a = CommandLine.arguments[1]

if CommandLine.arguments.count >= 3 {
    outBaseName = CommandLine.arguments[2]
}
if CommandLine.arguments.count == 4 && CommandLine.arguments[3] == "verbose" {
    verbose = true
}

let parseur = XMLGraphParse(anURL: URL(fileURLWithPath: a), verbose: verbose)
parseur.apply()
var translate = Dictionary<String, Int>()


var vertices = parseur.myVertices
for v in 0..<vertices.count {
    translate[vertices[v].id!] = v
}
let url1 = URL( fileURLWithPath: "\(outBaseName)\(outBaseName=="" ? "" : "_")vertex.dat" )
var content = ""
for v in vertices {
    if v.z != nil {
        content += "\(v.x!) \(v.y!) \(v.z!) \n"
    }else{
        content += "\(v.x!) \(v.y!) \n"
    }
}
try! content.write(to: url1, atomically: true, encoding: .utf8)


var edges = parseur.myEdges
let url2 = URL( fileURLWithPath: "\(outBaseName)\(outBaseName=="" ? "" : "_")edges.dat" )
content = ""
for e in edges {
    content += "\(translate[e.idA!]!) \(translate[e.idb!]!) \(e.flow==nil ? "" : "\(e.flow!)") \n"
}
try! content.write(to: url2, atomically: true, encoding: .utf8)


// export radius for each vertex
let url3 = URL( fileURLWithPath: "\(outBaseName)\(outBaseName=="" ? "" : "_")radius.dat" )
content = ""
for e in edges {
    if e.radius! > vertices[translate[e.idA!]!].radius ?? 0 {
        vertices[translate[e.idA!]!].radius = e.radius!
    }
    if e.radius! > vertices[translate[e.idb!]!].radius ?? 0 {
        vertices[translate[e.idb!]!].radius = e.radius!
    }
}
for v in 0..<vertices.count {
    content += "\(vertices[v].radius!)\n"
}

try! content.write(to: url3, atomically: true, encoding: .utf8)

print("converting done, graph exported in \(url1.lastPathComponent) \(url2.lastPathComponent) \(url3.lastPathComponent)")
