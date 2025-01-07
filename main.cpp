#include "Simulation.h"
#include "Vertex.h"
#include "Edge.h"
#include "Polygon.h"
#include "Cell.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
void readVTKAndCreateObjects(const std::string& filePath, 
                             std::vector<Vertex*>& vertices, 
                             std::vector<Edge*>& edges, 
                             std::vector<Polygon*>& polygons) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filePath);
    }

    std::string line;
    bool readingPoints = false;
    bool readingPolygons = false;
    std::vector<std::array<double, 3>> points;
    std::vector<std::vector<int>> polygonFaces;
    std::string keyword;
    int numpoints = -1;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!readingPoints && !readingPolygons){
            iss >> keyword;
        }
        if (keyword == "POINTS") {
            iss >> numpoints;
            readingPoints = true;
            readingPolygons = false;
            keyword = "reading points...";
            continue;
        } else if (keyword == "POLYGONS") {
            readingPoints = false;
            readingPolygons = true;
            keyword = "reading polygons...";
            continue;
        }

        if (readingPoints) {
            std::array<double, 3> point;
            iss >> point[0]>> point[1] >> point[2];
            points.push_back(point);
            numpoints--;
            if(numpoints==0){
                readingPoints = false;
            }
        } else if (readingPolygons) {
            int numVertices;
            iss >> numVertices;
            std::vector<int> face;
            for (int i = 0; i < numVertices; ++i) {
                int vertexId;
                iss >> vertexId;
                face.push_back(vertexId);
            }
            polygonFaces.push_back(face);
        }
    }
    file.close();

    // Create Vertex objects
    for (size_t i = 0; i < points.size(); ++i) {
        vertices.push_back(new Vertex(points[i], i));
    }

    // Create Edge objects
    std::map<std::tuple<int, int>, Edge*> edgeMap;
    int edgeId = 0;

    for (const auto& face : polygonFaces) {
        for (size_t i = 0; i < face.size(); ++i) {
            int v1 = face[i];
            int v2 = face[(i + 1) % face.size()];

            std::tuple edgeKey = std::make_tuple(v1, v2);
            std::tuple other = std::make_tuple(v2,v1);
            if (edgeMap.find(edgeKey) == edgeMap.end() && edgeMap.find(other)==edgeMap.end()) {
                edgeMap[edgeKey] = new Edge({vertices[v1], vertices[v2]}, edgeId++);
                edges.push_back(edgeMap[edgeKey]);
            }
        }
    }

    // Create Polygon objects
    int polygonId = 0;
    for (const auto& face : polygonFaces) {
        std::vector<Vertex*> faceVertices;
        std::vector<Edge*> faceEdges;
        for (size_t i = 0; i < face.size(); ++i) {
            faceVertices.push_back(vertices[face[i]]);

            int v1 = face[i];
            int v2 = face[(i + 1) % face.size()];

            std::tuple edgeKey = std::make_tuple(v1, v2);
            std::tuple other = std::make_tuple(v2, v1);
            if (edgeMap.find(edgeKey)!=edgeMap.end()){
                faceEdges.push_back(edgeMap[edgeKey]);
            }
            else if (edgeMap.find(other)!=edgeMap.end()){
                faceEdges.push_back(edgeMap[other]);
            }
            else {
                std::cout<<"Failed to find edge";
            }
        }

        polygons.push_back(new Polygon(faceVertices, faceEdges, polygonId++));        
    }
}


int main() {
    std::vector<Vertex*> vertices;
    std::vector<Edge*> edges;
    std::vector<Polygon*> polygons;


    // Step 4: Create the cell
    std::string shape = "cube";
    double eta = 1.0;
    double timestep = 0.001; //HY uses 0.005
    int numTimesteps = 5; //HY using 1000
    double V0 = 2.744; //HY using cube of length 1.4
    double A0 = 11.76;
    int log = 1; //HY using 10
    bool write = false;

    std::string vtkFilePath = std::string("vtk_in//")+shape+std::string(".vtk"); 
    readVTKAndCreateObjects(vtkFilePath, vertices, edges, polygons);
    std::vector<Cell*> cells;
    cells.push_back(new Cell(vertices, polygons, 0, V0, A0));
    std::cout << "Starting Simulation\n"; 
    Simulation sim(cells, polygons, edges, vertices, timestep, numTimesteps, eta, log, write); 

    // Cleanup dynamically allocated objects
    for (auto vertex : vertices) delete vertex;
    for (auto edge : edges) delete edge;
    for (auto polygon : polygons) delete polygon;
    for (auto cell:cells) delete cell;

    return 0;
}