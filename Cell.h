#pragma once
#include "Polygon.h" 
#include "Edge.h" 
#include "Vertex.h"
class Cell {
private:
    int id_;
    double volume_;
    double area_;
    std::array<double, 3> centroid_;
    std::vector<Polygon*> polygons_;
    std::vector<Vertex*> vertices_; 
    void updateVolume();
    void updateArea();
    void updateCentroid();
public:
    explicit Cell(std::vector<Vertex*> vertices, std::vector<Polygon*> polygons, int id);
    const double getVolume() const;
    const double getArea() const; 
    const std::array<double, 3>& getCentroid() const;
    const int getId() const;
    const std::vector<Polygon*>& getPolygons() const;
    const std::vector<Vertex*>& getVertices() const;
    void update();
    void setId(int id);

};
