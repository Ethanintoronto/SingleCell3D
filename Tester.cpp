#include <iostream>
#include <vector>
#include "Simulation.h"
#include "Vertex.h"
#include "Edge.h"
#include "Polygon.h"
#include "Cell.h"

class Tester {
public:
    // Constructor
    Tester() {}

    // Run all tests
    void runTests() {
        //testVertex();
        //testEdge();
        testPolygon();
        //testCell();
        //testSimulation();
    }

private:
    // Test Vertex functionality
    void testVertex() {
        std::cout << "Running Vertex tests...\n";

        // Create a Vertex and verify its properties
        std::array<double, 3> position = {1.0, 2.0, 3.0};
        Vertex v(position, 0);

        // Check position values
        const auto& pos = v.getPos();
        if (pos[0] == 1.0 && pos[1] == 2.0 && pos[2] == 3.0) {
            std::cout << "Vertex position test passed.\n";
        } else {
            std::cerr << "Vertex position test failed.\n";
        }

        // Check ID
        if (v.getId() == 0) {
            std::cout << "Vertex ID test passed.\n";
        } else {
            std::cerr << "Vertex ID test failed.\n";
        }
    }

    // Test Edge functionality
    void testEdge() {
        std::cout << "Running Edge tests...\n";

        // Create vertices for testing edges
        std::array<double, 3> pos1 = {0.0, 0.0, 0.0};
        std::array<double, 3> pos2 = {1.0, 0.0, 0.0};
        Vertex v1(pos1, 1);
        Vertex v2(pos2, 2);

        // Create an Edge and verify its properties
        std::array<Vertex*, 2> vertices = {&v1, &v2};
        //Edge e(vertices, 0);
        Edge * e = new Edge(vertices, 0);
        // Check the position of the vertices in the Edge:
        if (e->getVertices()[0]->getPos()[0] ==0 &&e->getVertices()[1]->getPos()[0] ==1 ){
            std::cout << "Edge position test passed.\n";
        }
        else {
            std::cerr << "Edge position test failed. \n";
        }
        // Check connected vertices
        if (e->getVertices()[0]->getId() == 1 && e->getVertices()[1]->getId() == 2) {
            std::cout << "Edge connectivity test passed.\n";
        } else {
            std::cerr << "Edge connectivity test failed.\n";
        }
        // Check edge update function
        e->update();
        if (e->getCentroid()[0]==0.5){
            std::cout<<"centroid test passed.\n";
        } 
        else {
            std::cout<<"Centroid: "<<e->getCentroid()[0]<<","<<e->getCentroid()[1]<<","<<e->getCentroid()[2]<<"\n"; 
            std::cout<<"Centroid test failed\n";
        }
    }

    // Test Polygon functionality
    void testPolygon() {
    std::cout << "Running Polygon tests...\n";

    // Cube vertices: a cube with side length 1 centered at the origin
    std::array<std::array<double, 3>, 8> cubeVertices = {{
        {{-0.5, -0.5, -0.5}}, {{ 0.5, -0.5, -0.5}}, {{ 0.5,  0.5, -0.5}}, {{-0.5,  0.5, -0.5}},
        {{-0.5, -0.5,  0.5}}, {{ 0.5, -0.5,  0.5}}, {{ 0.5,  0.5,  0.5}}, {{-0.5,  0.5,  0.5}}
    }};
    
    // Step 1: Create the vertices for the cube
    std::vector<Vertex*> vertices;
    for (int i = 0; i < 8; ++i) {
        vertices.push_back(new Vertex(cubeVertices[i], i));
    }

    // Step 2: Define the edges (each edge connects two vertices)
    std::vector<Edge*> edges;

    // Define cube edges using pairs of vertices (e.g., edge from vertex 0 to vertex 1)
    std::array<std::array<int, 2>, 12> edgeIndices = {{
        {{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}}, // Bottom face edges
        {{4, 5}}, {{5, 6}}, {{6, 7}}, {{7, 4}}, // Top face edges
        {{0, 4}}, {{1, 5}}, {{2, 6}}, {{3, 7}}  // Vertical edges
    }};

    for (int i = 0; i < 12; ++i) {
        std::array<Vertex*, 2> edgeVertices = {vertices[edgeIndices[i][0]], vertices[edgeIndices[i][1]]};
        edges.push_back(new Edge(edgeVertices, i));
    }
    // Step 3: Define the polygons (cube has 6 faces, each face is a square)
    std::vector<Polygon*> polygons;

    // Define each face as a polygon with its vertices and edges
    std::array<std::array<int, 4>, 6> polygonVIndices = {{
        {{0, 1, 2, 3}}, {{4, 5, 6, 7}}, // Bottom and Top faces
        {{0, 1, 5, 4}}, {{1, 2, 6, 5}}, // Side faces
        {{2, 3, 7, 6}}, {{3, 0, 4, 7}}  // Side faces
    }};
    std::array<std::array<int, 4>, 6> polygonEIndices = {{
        {{0, 1, 2, 3}}, {{4, 5, 6, 7}}, // Bottom and Top faces
        {{0, 9, 4, 8}}, {{1, 10, 5, 9}}, // Side faces
        {{2, 11, 6, 10}}, {{3, 8, 7, 11}}  // Side faces
    }};

    for (int i = 0; i < 6; ++i) {
        std::vector<Vertex*> polygonVertices;
        std::vector<Edge*> polygonEdges;
        
        // Get the vertices and edges for each polygon
        for (int j = 0; j < 4; ++j) {
            polygonVertices.push_back(vertices[polygonVIndices[i][j]]);
        }
        // Add edges corresponding to each face (assuming that each face is a quad)
        for (int j = 0; j < 4; ++j) {
            polygonEdges.push_back(edges[polygonEIndices[i][j]]);
        }
        polygons.push_back(new Polygon(polygonVertices, polygonEdges, i));
    }

        // Check polygon vertex count
        if (polygons[0]->getVertices().size() == 4) {
            std::cout << "Polygon vertex count test passed.\n";
        } else {
            std::cerr << "Polygon vertex count test failed.\n";
        }
        //Test Polygon update function
        //polygons[0]->update();
        std::cout <<"Centroid: (" << polygons[0]->getCentroid()[0]<<","<<polygons[0]->getCentroid()[1] <<", "<<polygons[0]->getCentroid()[2]<< ")\n";
        std::cout <<"Area: " <<polygons[0]->getArea()<<"\n";
        std::cout <<"Area Vector: ("<<polygons[0]->getAreaVector()[0]<<","<<polygons[0]->getAreaVector()[1]<<","<<polygons[0]->getAreaVector()[2]<<")\n"; 
        std::cout <<"Perimeter: " <<polygons[0]->getPerimeter()<<"\n";
        std::cout <<"Edge centroid: ("<<polygons[0]->getEdges()[2]->getCentroid()[0]<<","<<polygons[0]->getEdges()[2]->getCentroid()[1]<<","<<polygons[0]->getEdges()[2]->getCentroid()[2]<<")\n";
        std::cout <<"Edge length: "<<polygons[0]->getEdges()[3]->getLength()<<"\n";
    }

    // Test Cell functionality
    void testCell() {
        std::cout << "Running Cell tests...\n";

        // Create vertices, edges, and a polygon to form a cell
        std::array<double, 3> pos1 = {0.0, 0.0, 0.0};
        std::array<double, 3> pos2 = {1.0, 0.0, 0.0};
        std::array<double, 3> pos3 = {1.0, 1.0, 0.0};
        std::array<double, 3> pos4 = {0.0, 1.0, 0.0};

        Vertex v1(pos1, 1), v2(pos2, 2), v3(pos3, 3), v4(pos4, 4);
        std::vector<Vertex*> vertices = {&v1, &v2, &v3, &v4};

        std::array<Vertex*, 2> edgeVertices1 = {&v1, &v2};
        std::array<Vertex*, 2> edgeVertices2 = {&v2, &v3};
        std::array<Vertex*, 2> edgeVertices3 = {&v3, &v4};
        std::array<Vertex*, 2> edgeVertices4 = {&v4, &v1};

        Edge e1(edgeVertices1, 1), e2(edgeVertices2, 2), e3(edgeVertices3, 3), e4(edgeVertices4, 4);
        std::vector<Edge*> edges = {&e1, &e2, &e3, &e4};

        Polygon p(vertices, edges, 0);
        std::vector<Polygon*> polygons = {&p};

        // Create a Cell
        Cell * c = new Cell(vertices, polygons, 1);

        // Check properties of the cell
        if (c->getId() == 1) {
            std::cout << "Cell ID test passed.\n";
        } else {
            std::cerr << "Cell ID test failed.\n";
        }

        //Check Cell.update()
        c->update();

    }

    // Test Simulation functionality
    void testSimulation() {
        std::cout << "Running Simulation tests...\n";

        // Create vertices, edges, and a polygon to form a cell
        std::array<double, 3> pos1 = {0.0, 0.0, 0.0};
        std::array<double, 3> pos2 = {1.0, 0.0, 0.0};
        std::array<double, 3> pos3 = {1.0, 1.0, 0.0};
        std::array<double, 3> pos4 = {0.0, 1.0, 0.0};

        Vertex v1(pos1, 1), v2(pos2, 2), v3(pos3, 3), v4(pos4, 4);
        std::vector<Vertex*> vertices = {&v1, &v2, &v3, &v4};

        std::array<Vertex*, 2> edgeVertices1 = {&v1, &v2};
        std::array<Vertex*, 2> edgeVertices2 = {&v2, &v3};
        std::array<Vertex*, 2> edgeVertices3 = {&v3, &v4};
        std::array<Vertex*, 2> edgeVertices4 = {&v4, &v1};

        Edge e1(edgeVertices1, 1), e2(edgeVertices2, 2), e3(edgeVertices3, 3), e4(edgeVertices4, 4);
        std::vector<Edge*> edges = {&e1, &e2, &e3, &e4};

        Polygon p(vertices, edges, 0);
        std::vector<Polygon*> polygons = {&p};

        Cell c(vertices, polygons, 1);
        std::vector<Cell*> cells = {&c};

        double timestep = 0.01;
        int numTimesteps = 1;
        double Kv = 1.0, Ka = 1.0, V0 = 1.0, A0 = 1.0, eta = 1.0;
        int log = 1;
        bool write = true;

        Simulation sim(cells, polygons, edges, vertices, timestep, numTimesteps, Kv, Ka, V0, A0, eta, log, write);

        // Check the results (adjust based on your simulation logic)
        std::cout << "Simulation test completed.\n";

        // Clean up allocated memory
        for (auto c : cells) delete c;
    }
};

int main2() {
    Tester tester;
    tester.runTests();
    return 0;
}
