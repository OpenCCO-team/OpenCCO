//
//  main.swift
//  Graph2statBifRad
//
//  Created by Bertrand Kerautret on 30/01/2023.
//

import Foundation

struct Vertex {
    var id: Int? = nil
    var x: Double? = nil
    var y: Double? = nil
    var z: Double? = nil
    var radius: Double? = nil
    var children1: Int? = nil
    var children2: Int? = nil
    var parent: Int? = nil
    var bifLevel: Int = 0
}

struct Edge {
  var idA: Int?
  var idB: Int?
  var radius: Double?
}



func readEdges(fileName: String) ->  [Edge]
{
    var result: [Edge] = []
    if let path = URL(string: fileName)?.absoluteString {
        do {
            let data = try String(contentsOfFile: path, encoding: .utf8)
            let myStrings = data.components(separatedBy: .newlines)
            for l in myStrings {
                let elem = l.components(separatedBy: " ")
                if elem.count > 1 {
                    if let ind1 = Int(elem[0]), let ind2 = Int(elem[1]){
                        let e = Edge(idA: ind1, idB: ind2)
                        result.append(e)
                    }
                }
            }
        } catch {
            print(error)
        }
    }
    return result
}



func readVertex(fileName: String) ->  [Vertex]
{
    var result: [Vertex] = []
    var nb = 0
    if let path = URL(string: fileName)?.absoluteString {
        do {
            let data = try String(contentsOfFile: path, encoding: .utf8)
            let myStrings = data.components(separatedBy: .newlines)
            for l in myStrings {
                let elem = l.components(separatedBy: " ")
                if elem.count > 2  {
                    if let x = Double (elem[0]),
                        let y = Double(elem[1]),
                        let z = Double(elem[2])
                    {
                        result.append(Vertex(id: nb, x: x, y: y, z:z))
                        nb += 1
                    }
                }
            }
        } catch {
            print(error)
        }
    }
    return result
}



func readRadius(fileName: String) ->  [Double]
{
    var result: [Double] = []
    if let path = URL(string: fileName)?.absoluteString {
        do {
            let data = try String(contentsOfFile: path, encoding: .utf8)
            let myStrings = data.components(separatedBy: .newlines)
            for l in myStrings {
                let elem = l.components(separatedBy: " ")
                if elem.count > 0  {
                    if let r = Double (elem[0])
                    {
                        result.append(r)
                    }
                }
            }
        } catch {
            print(error)
        }
    }
    return result
}


func getMeanVariance(t: [Double]) -> (Double, Double) {
    var mean: Double = 0.0
    var variance: Double = 0.0
    for v in t {
        mean += v
    }
    mean /= (Double)(t.count)
    for v in t {
        variance += (v-mean)*(v-mean)
    }
    variance /= (Double)(t.count)
    return (mean, variance)
}

func usage(){
    print("Command line args: \(CommandLine.arguments[0]) nodesFile edgeFile radiusStatFile resultFile")
}


if CommandLine.arguments.count < 5 {
    usage()
    exit(1)
}

let nodesFile = CommandLine.arguments[1]
let edgeFile = CommandLine.arguments[2]
let radiusFile = CommandLine.arguments[3]
let radiusStatFile = CommandLine.arguments[4]

print("Reading edges files \(edgeFile)", terminator: "")
var edges = readEdges(fileName: edgeFile)
print("[done]: \(edges.count)")


print("Reading nodes files \(nodesFile)", terminator: "")
var vertex = readVertex(fileName: nodesFile)
print("[done]: \(vertex.count)")


print("Reading radius files \(radiusFile)", terminator: "")
let radius = readRadius(fileName: radiusFile)
print("[done]: \(radius.count)")



// updating for each vertex, the radius
for v in vertex {
    if let i = v.id {
        vertex[i].radius = radius[i]
        vertex[i].bifLevel = -1
    }
}


// Updating for each vertex, the childs 1 and 2 from the edges
for e in edges {
    if let vA = e.idA, let vB = e.idB {
        if vertex[vA].children1 == nil {
            vertex[vA].children1 = vB
        }else {
            vertex[vA].children2 = vB
        }
        vertex[vB].parent =  vA
    }
}

// Update for each vertex its bifurcation level
var listTraite: [Int] = [0]
while !listTraite.isEmpty
{
    let s = listTraite.removeLast()
    if let p = vertex[s].parent  {
        vertex[s].bifLevel = vertex[p].bifLevel + 1
    }else {
       vertex[s].bifLevel = 0
    }
    if let c = vertex[s].children1 {
        listTraite.append(c)
    }
    if let c = vertex[s].children2 {
        listTraite.append(c)
    }
}



var bifLevelMax = 0
for v in vertex {
    if v.bifLevel > bifLevelMax {
        bifLevelMax = v.bifLevel
    }
}


// Initialisation of the tab representing for each bifurcation level le tabular of radius of the edges.
var tabRadiusBif = [[Double]]()
for _ in 0...bifLevelMax {
    tabRadiusBif.append([Double]())
}

// Filling tab radius
for v in vertex {
    if let r = v.radius {
        tabRadiusBif[v.bifLevel].append(r)
    }
}

let url = URL( fileURLWithPath: radiusStatFile )
var content = ""
content += "# Mean Variance nbElements \n"

// Starting from level 1 since the root is by definition not associated to any segment
for i in 1...bifLevelMax {
 let meanEsp = getMeanVariance(t: tabRadiusBif[i])
    content += "\(i-1) \(meanEsp.0) \(meanEsp.1)  \(tabRadiusBif[i].count) \n"
}

try! content.write(to: url, atomically: true, encoding: .utf8)
