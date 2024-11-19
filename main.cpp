#include "Simulation.h"
#include "Vertex.h"
#include "Edge.h"
#include "Polygon.h"
#include "Cell.h"
int main() {
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
    std::array<std::array<int, 4>, 6> polygonIndices = {{
        {{0, 1, 2, 3}}, {{4, 5, 6, 7}}, // Bottom and Top faces
        {{0, 1, 5, 4}}, {{1, 2, 6, 5}}, // Side faces
        {{2, 3, 7, 6}}, {{3, 0, 4, 7}}  // Side faces
    }};

    for (int i = 0; i < 6; ++i) {
        std::vector<Vertex*> polygonVertices;
        std::vector<Edge*> polygonEdges;
        
        // Get the vertices and edges for each polygon
        for (int j = 0; j < 4; ++j) {
            polygonVertices.push_back(vertices[polygonIndices[i][j]]);
        }
        // Add edges corresponding to each face (assuming that each face is a quad)
        for (int j = 0; j < 4; ++j) {
            int idx1 = polygonIndices[i][j];
            int idx2 = polygonIndices[i][(j + 1) % 4];
            polygonEdges.push_back(edges[idx1 * 2 + idx2 % 2]);  // This line assumes edge indexing
        }
        polygons.push_back(new Polygon(polygonVertices, polygonEdges, i));
    }

    // Step 4: Create the cell
    std::vector<Cell*> cells;
    cells.push_back(new Cell(vertices, polygons, 1));
    double eta = 1.0;
    double timestep = 0.01;
    int numTimesteps = 10;
    double Kv = 10.0;
    double Ka = 1.0;
    double V0 = 1.0;
    double A0 = 1.0;
    int log = 1; 

    Simulation sim(cells, vertices, timestep, numTimesteps, Kv, Ka, V0, A0, eta, log); 

    // Cleanup dynamically allocated memory
    for (auto v : vertices) delete v;
    for (auto e : edges) delete e;
    for (auto p : polygons) delete p;

    return 0;
}