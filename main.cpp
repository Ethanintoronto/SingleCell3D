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
        Polygon* polygon = new Polygon(faceVertices, faceEdges, polygonId++);
        polygons.push_back(polygon);
        //polygons.push_back(new Polygon(faceVertices, faceEdges, polygonId++));        
    }
}

int main() {
    std::vector<Vertex*> vertices;
    std::vector<Edge*> edges;
    std::vector<Polygon*> polygons;

    // Initialize Simulation Parameters:
    std::string shape = "cube";
    double eta = 1.0;
    double timestep = 0.001; 
    int numTimesteps = 10000; 
    double V0 = 1; 
    double A0 = 6;
    double gamma = 5;
    double Kv = 10.;
    double KaCell = 0;
    double KaPoly = 10.;

    int log = 10; 
    bool write = true;

    std::string vtkFilePath = std::string("vtk_in//")+shape+std::string(".vtk"); 
    readVTKAndCreateObjects(vtkFilePath, vertices, edges, polygons);
    
    //Set gamma parameters
    for (Polygon* polygon :polygons){
        polygon->setKa(KaPoly);
        if (polygon->getId() == 2){
            polygon->setGamma(gamma);
        }
        else if (polygon->getId() == 5){
            polygon->setGamma(-1*gamma);
        }
    }
    std::vector<Cell*> cells;
    Cell* cell = new Cell(vertices, polygons, 0, V0, A0);
    cell->setKv(Kv);
    cell->setKa(KaCell); 
    cells.push_back(cell);
    std::cout << "Starting Simulation\n"; 
    Simulation sim(cells, polygons, edges, vertices, timestep, numTimesteps, eta, log, write); 

    // Cleanup dynamically allocated objects
    for (auto vertex : vertices) delete vertex;
    for (auto edge : edges) delete edge;
    for (auto polygon : polygons) delete polygon;
    for (auto cell:cells) delete cell;

    return 0;
}