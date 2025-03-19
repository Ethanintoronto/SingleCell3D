#include "Simulation.h"
#include "Vertex.h"
#include "Edge.h"
#include "Polygon.h"
#include "Cell.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <filesystem>
#include <regex>
#include <ctime>
#include <unordered_set>
void readVTKAndCreateObjects(const std::string& filePath, 
                             std::vector<Vertex*>& vertices, 
                             std::vector<Edge*>& edges, 
                             std::vector<Polygon*>& polygons,
                             std::vector<Cell*>& cells,
                             double V0, double A0) { // Group polygons by cube_id
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
    int numpoints = -1, numPolygons = -1, numPolygonsLeft = -1;

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
            numPolygonsLeft = numPolygons;
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
            numPolygonsLeft--;
            if (numPolygonsLeft == 0){
                readingPolygons = false;
            }
        }
    }
    file.close();

    //Open helper to read cells
    std::ifstream helperFile(filePath.substr(0,filePath.length()-4)+"_helper.txt");
    if (!helperFile.is_open()) {
        throw std::runtime_error("Could not open file " + filePath.substr(0,filePath.length()-4)+"_helper.txt");
    }
    bool readingCells = false;
    int polyId;
    std::vector<std::vector<int>> cellStructure;
    while (std::getline(helperFile, line)) {
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
    helperFile.close();
    

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

    //Assign boundary designation to boundary polygons:
    std::unordered_set<int> seenPolyIds;
    std::vector<bool> boundPolys(numPolygons, true);
    for (const auto& cell: cellStructure){  
        for (int poly_id: cell){
            if (seenPolyIds.count(poly_id)>0){
                boundPolys[poly_id] = false;
            }
            else{
                seenPolyIds.insert(poly_id);
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
        Polygon* polygon = new Polygon(faceVertices, faceEdges, polygonId);
        if (boundPolys[polygonId]){
            polygon->setBoundary(true);
        }
        polygons.push_back(polygon);  
        polygonId++;     
    }

    // Create cells
    int cellId = 0;
    for (int i = 0; i<cellStructure.size();i++){
        const auto& thisCellStructure = cellStructure[i];
        std::vector<Vertex*> cellVertices;
        std::vector<Polygon*> cellPolygons;
        for (int j = 0; j<thisCellStructure.size(); j++){
            Polygon* polygon = polygons[thisCellStructure[j]];
            cellPolygons.push_back(polygon);
            for (int k = 0; k<polygon->getVertices().size();k++){
                Vertex* vertex = polygon->getVertices()[k];
                if (std::find(cellVertices.begin(), cellVertices.end(), vertex) == cellVertices.end()) { //Vertex not already added
                    cellVertices.push_back(vertex);
                }
            }
        }
        Cell* cell = new Cell(cellVertices, cellPolygons, cellId++, V0, A0);
        cells.push_back(cell);
    }
    
}
int get_next_index(const std::string& directory) {
    std::regex pattern(R"(_(\d{3})$)"); // Matches "_000", "_001", etc. at the end of filenames
    if (!std::filesystem::exists(directory) || std::filesystem::is_empty(directory)) {
        return 0; // Return 0 if the directory does not exist or is empty
    }
    int max_index = -1;

    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_directory()) {
            std::smatch match;
            std::string dirname = entry.path().stem().string(); // Get dirname without extension

            if (std::regex_search(dirname, match, pattern)) {
                int num = std::stoi(match[1]); // Extract integer
                max_index = std::max(max_index, num);
            }
        }
    }

    return max_index + 1;
}

std::string getDate() {
    // Get the current time
    std::time_t t = std::time(nullptr);

    // Convert it to a tm structure
    std::tm* now = std::localtime(&t);

    // Create a string stream to format the date
    std::ostringstream date_stream;
    date_stream << std::put_time(now, "%Y-%m-%d"); // Format: YYYY-MM-DD

    // Get the formatted date as a string
    std::string date_str = date_stream.str();

    return date_str;
}


int main() {
    bool batchMode = false;
    if (batchMode){
        int period_inc = 50; // 1 tau = 0.5 min, inc +50 = +0.5 tau = + 0.25 min  
        double gamma = 0.5;
        double gamma_inc = 0.25;
        int num_periods = 10;
        int num_gammas = 10;

        double Ks = 0.5;
        double Ks_inc = 0.5;
        double KsTrailing_inc = 0.5;
        for (int i = 0; i<num_gammas; i++){
            
            
            //int period = 250;
            double KsTrailing = Ks;


            for (int j = 0; j<num_periods; j++){
                std::vector<Vertex*> vertices;
                std::vector<Edge*> edges;
                std::vector<Polygon*> polygons;
                std::vector<Cell*> cells;
                // Initialize Simulation Parameters:
                std::string shape = "stacked_cubes";
                std::string directory = "data/" + getDate(); 
                int id = get_next_index(directory);
                double mu = 1.0;
                double V0 = 1.; 
                double A0 = 6.;
                double Kv = 10.;
                double Ka = 1.;
                //double Ks = 1.;
                double gamma = 2;
                int period = 500;
                //double KsTrailing = 10;
                bool boundary = false;
                int midSteps = 0;

                //period = 500 timesteps = 5 tau = 2.5 min
                double tau = 1/(mu*Kv*V0); //0.1 
                double timestep = 0.01*tau; 

                //run for 4 periods + 60 tau 
                //int numTimesteps = period*4 + 100*60; //100 integration timesteps/tau * 60 tau
                int numTimesteps = 2000 + 100*60;
                int log = period/10; 
                bool write = true;

                std::string vtkFilePath = std::string("vtk_in//")+shape+std::string(".vtk"); 
                readVTKAndCreateObjects(vtkFilePath, vertices, edges, polygons, cells, V0, A0);
                //Set gamma parameters
                for (Polygon* polygon :polygons){

                    polygon->setKs(Ks);
                    if (polygon->getId() == 1 ||polygon->getId() == 7||polygon->getId() == 12||polygon->getId() == 17){
                        polygon->setGamma(gamma);
                    }
                    else if (polygon->getId() == 3||polygon->getId() == 9 ||polygon->getId() == 14||polygon->getId() == 19){
                        polygon->setGamma(-1*gamma);
                        polygon->setKs(KsTrailing);
                    }
                }
                for (Cell* cell : cells){
                    cell->setKv(Kv);
                    cell->setKa(Ka);
                }
                std::cout << "Starting Simulation: " <<std::setfill('0') <<std::setw(3) << id<<std::endl; 
                Simulation sim(cells, polygons, edges, vertices, id, period, timestep, numTimesteps, mu, log, write); 
                sim.setBoundary(boundary);
                sim.setMidSteps(midSteps);
                sim.Run();

                // Cleanup dynamically allocated objects
                for (auto vertex : vertices) delete vertex;
                for (auto edge : edges) delete edge;
                for (auto polygon : polygons) delete polygon;
                for (auto cell:cells) delete cell;


                //period += period_inc;
                KsTrailing+=KsTrailing_inc;
            }


            //gamma+= gamma_inc;
            Ks+= Ks_inc;
        }
    }
    else if (!batchMode){
        std::vector<Vertex*> vertices;
        std::vector<Edge*> edges;
        std::vector<Polygon*> polygons;
        std::vector<Cell*> cells;
        // Initialize Simulation Parameters:

        std::string shape = "stacked_cubes";
        std::string directory = "data/" + getDate(); 
        int id = get_next_index(directory);
        double mu = 1.0;
        double V0 = 1; 
        double A0 = 6;
        double gamma = 2.;
        double Kv = 10.;
        double Ka = 1.;
        double Ks = 1.;
        double KsTrailing = 1.;
        double period = 500;
        bool boundary = false;
        int midSteps = 0;

        double tau = 1/(mu*Kv*V0); //0.1 
        double timestep = 0.01*tau;

        //run for 4 periods + 60 tau 
        int numTimesteps = period*4 + 100*60; //timesteps/log = timesteps/period * periods/log = 500 * 1/10 

        int log = period/10; 
        bool write = true;

        std::string vtkFilePath = std::string("vtk_in//")+shape+std::string(".vtk"); 
        readVTKAndCreateObjects(vtkFilePath, vertices, edges, polygons,cells, V0, A0);
        
        //Set gamma parameters
        for (Polygon* polygon :polygons){
            polygon->setKs(Ks);
            if (polygon->getId() == 1 ||polygon->getId() == 7||polygon->getId() == 12||polygon->getId() == 17){
                polygon->setGamma(gamma);
            }
            else if (polygon->getId() == 3||polygon->getId() == 9 ||polygon->getId() == 14||polygon->getId() == 19){
                polygon->setGamma(-1*gamma);
                polygon->setKs(KsTrailing);
            }
        }
        for (Cell* cell : cells){
            cell->setKv(Kv);
            cell->setKa(Ka);
        }
        std::cout << "Starting Simulation: " <<std::setfill('0') <<std::setw(3) << id<<std::endl; 
        Simulation sim(cells, polygons, edges, vertices, id, period, timestep, numTimesteps, mu, log, write); 
        sim.setBoundary(boundary);
        sim.setMidSteps(midSteps);
        sim.Run();
        // Cleanup dynamically allocated objects
        for (auto vertex : vertices) delete vertex;
        for (auto edge : edges) delete edge;
        for (auto polygon : polygons) delete polygon;
        for (auto cell: cells) delete cell;
    }
    return 0;
}