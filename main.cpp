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
    }
}

int main() {
    int id = 106;
    bool batchMode = false;
    if (batchMode){
        int id = 6;
        double period_inc = 50; 
        double gamma = 1;
        double gamma_inc = 1;
        for (int i = 0; i<10; i++){
            double period = 250;
            for (int j = 0; j<10; j++){
                std::vector<Vertex*> vertices;
                std::vector<Edge*> edges;
                std::vector<Polygon*> polygons;
                // Initialize Simulation Parameters:
                std::string shape = "cube";
                double mu = 1.0;
                double V0 = 1.; 
                double A0 = 6.;
                double Kv = 10.;
                double Ka = 1.;
                double Ks = 1.;

                //period = 500 timesteps = 5 tau = 2.5 min
                double tau = 1/(mu*Kv*V0); //0.1 
                double timestep = 0.01*tau; 

                //run for 4 periods + 60 tau 
                //int numTimesteps = period*4 + 100*60; //100 integration timesteps/tau * 60 tau
                int numTimesteps = 2000 + 100*60;
                int log = period/10; 
                bool write = true;

                std::string vtkFilePath = std::string("vtk_in//")+shape+std::string(".vtk"); 
                readVTKAndCreateObjects(vtkFilePath, vertices, edges, polygons);
                
                //Set gamma parameters
                for (Polygon* polygon :polygons){
                    polygon->setKs(Ks);
                    if (polygon->getId() == 2){
                        polygon->setGamma(gamma);
                    }
                    else if (polygon->getId() == 5){
                        polygon->setGamma(-1*gamma);
                        polygon->setKs(10.);
                    }
                }
                std::vector<Cell*> cells;
                Cell* cell = new Cell(vertices, polygons, 0, V0, A0);
                cell->setKv(Kv);
                cell->setKa(Ka); 
                cells.push_back(cell);
                std::cout << "Starting Simulation\n"; 
                Simulation sim(cells, polygons, edges, vertices, id, period, timestep, numTimesteps, mu, log, write); 

                // Cleanup dynamically allocated objects
                for (auto vertex : vertices) delete vertex;
                for (auto edge : edges) delete edge;
                for (auto polygon : polygons) delete polygon;
                for (auto cell:cells) delete cell;

                id++;
                period += period_inc;
            }
            gamma+= gamma_inc;
        }
    }
    else{
        std::vector<Vertex*> vertices;
        std::vector<Edge*> edges;
        std::vector<Polygon*> polygons;
        // Initialize Simulation Parameters:
        std::string shape = "cube";
        double mu = 1.0;
        double V0 = 1.; 
        double A0 = 6.;
        double gamma = 2.;
        double Kv = 10.;
        double Ka = 1.;
        double Ks = 1.;
        double KsTrailing = 10.;
        double period = 250;

        //period = 500 timesteps = 5 tau = 2.5 min
        double tau = 1/(mu*Kv*V0); //0.1 
        double timestep = 0.01*tau; 

        //run for 4 periods + 60 tau 
        int numTimesteps = period*4 + 100*60; //100 integration timesteps/tau * 60 tau
        int log = period/10; 
        bool write = true;

        std::string vtkFilePath = std::string("vtk_in//")+shape+std::string(".vtk"); 
        readVTKAndCreateObjects(vtkFilePath, vertices, edges, polygons);
        
        //Set gamma parameters
        for (Polygon* polygon :polygons){
            polygon->setKs(Ks);
            if (polygon->getId() == 2){
                polygon->setGamma(gamma);
            }
            else if (polygon->getId() == 5){
                polygon->setGamma(-1*gamma);
                polygon->setKs(KsTrailing);
            }
        }
        std::vector<Cell*> cells;
        Cell* cell = new Cell(vertices, polygons, 0, V0, A0);
        cell->setKv(Kv);
        cell->setKa(Ka); 
        cells.push_back(cell);
        std::cout << "Starting Simulation\n"; 
        Simulation sim(cells, polygons, edges, vertices, id, period, timestep, numTimesteps, mu, log, write); 

        // Cleanup dynamically allocated objects
        for (auto vertex : vertices) delete vertex;
        for (auto edge : edges) delete edge;
        for (auto polygon : polygons) delete polygon;
        for (auto cell:cells) delete cell;
    }
    return 0;
}