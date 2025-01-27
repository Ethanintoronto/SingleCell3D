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
                             std::vector<Polygon*>& polygons) { // Group polygons by cube_id
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
    int numpoints = -1, numPolygons = -1;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!readingPoints && !readingPolygons) {
            iss >> keyword;
        }

        if (keyword == "POINTS") {
            iss >> numpoints;
            readingPoints = true;
            readingPolygons = false;
            keyword = "reading points...";
            continue;
        } else if (keyword == "POLYGONS") {
            iss >> numPolygons; // Read the number of polygons
            readingPoints = false;
            readingPolygons = true;
            keyword = "reading polygons...";
            continue;
        }

        if (readingPoints) {
            std::array<double, 3> point;
            iss >> point[0] >> point[1] >> point[2];
            points.push_back(point);
            numpoints--;
            if (numpoints == 0) {
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
            numPolygons--;
            if (numPolygons == 0){
                readingPolygons = false;
            }
        }
    }
    file.close();

    //Open helper to read cells
    std::ifstream file(filePath.substr(0,filePath.length()-4)+"_helper.txt");
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filePath);
    }
    std::string line;
    bool readingCells = false;
    int polyId;
    std::vector<std::vector<int>> cellStructure;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!readingCells) {
            readingCells = true;
            continue;
        }
        if (readingCells){
            iss >> numPolygons;
            std::vector<int> cell;
            for (int i = 0; i<numPolygons; i++){
                iss >> polyId;
                cell.push_back(polyId);
            }
            cellStructure.push_back(cell);
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
            std::tuple other = std::make_tuple(v2, v1);
            if (edgeMap.find(edgeKey) == edgeMap.end() && edgeMap.find(other) == edgeMap.end()) {
                edgeMap[edgeKey] = new Edge({vertices[v1], vertices[v2]}, edgeId++);
                edges.push_back(edgeMap[edgeKey]);
            }
        }
    }

    // Create Polygon objects and group them by cube_id
    int polygonId = 0;
    for (size_t i = 0; i < polygonFaces.size(); ++i) {
        const auto& face = polygonFaces[i];
        std::vector<Vertex*> faceVertices;
        std::vector<Edge*> faceEdges;

        for (size_t j = 0; j < face.size(); ++j) {
            faceVertices.push_back(vertices[face[j]]);

            int v1 = face[j];
            int v2 = face[(j + 1) % face.size()];
            
            std::tuple edgeKey = std::make_tuple(v1, v2);
            std::tuple other = std::make_tuple(v2, v1);
            if (edgeMap.find(edgeKey) != edgeMap.end()) {
                faceEdges.push_back(edgeMap[edgeKey]);
            } else if (edgeMap.find(other) != edgeMap.end()) {
                faceEdges.push_back(edgeMap[other]);
            } else {
                std::cerr << "Failed to find edge\n";
            }
        }
        // Create a polygon
        Polygon* polygon = new Polygon(faceVertices, faceEdges, polygonId++);
        polygons.push_back(polygon);
    }

    // Create cells
    for (int i = 0; i<cellStructure.size();i++){
        const auto& thisCellStructure = cellStructure[i];
        std::vector<Polygon*> cellPolygons;
        for (int j = 0; j<thisCellStructure.size(); j++){
            cellPolygons.push_back(polygons[thisCellStructure[j]]);
        }
        Cell* cell = new Cell(cellVertices, cellPolygons, cellId++, V0, A0);
    }
    
}


int main() {
    std::vector<Vertex*> vertices;
    std::vector<Edge*> edges;
    std::vector<Polygon*> polygons;
    std::vector<Cell*> cells;
    // Initialize Simulation Parameters:
    std::string file_in = "stacked_cubes.vtk";
    double eta = 1.0;
    double timestep = 0.001; 
    int numTimesteps = 10; 
    double V0 = 1; 
    double A0 = 6;
    double gamma = 0;
    double Kv = 10.;
    double KaCell = 0;
    double KaPoly = 10.;

    int log = 1; 
    bool write = true;

    std::string vtkFilePath = std::string("vtk_in//")+file_in; 
    readVTKAndCreateObjects(vtkFilePath, vertices, edges, polygons);
    
    //Build cubes

    std::cout << "Starting Simulation\n"; 
    Simulation sim(cells, polygons, edges, vertices, timestep, numTimesteps, eta, log, write); 

    // Cleanup dynamically allocated objects
    for (auto vertex : vertices) delete vertex;
    for (auto edge : edges) delete edge;
    for (auto polygon : polygons) delete polygon;
    for (auto cell:cells) delete cell;

    return 0;
}