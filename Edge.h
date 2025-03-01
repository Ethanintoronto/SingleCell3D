#pragma once
#include "Vertex.h"
class Edge{
private:
    int id_;
    double length_;
    std::array<double, 3> centroid_;
    std::array<Vertex*,2> vertices_;
    void updateLength();  
    void updateCentroid();
public:
    explicit Edge(std::array<Vertex*,2> vertices, int id);

    const std::array<Vertex*,2>& getVertices() const;

    int getId() const;

    double getLength() const; 

    const std::array<double, 3>& getCentroid() const;

    void update();

    void setId(int id);
};