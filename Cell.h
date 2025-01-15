#pragma once
#include "Polygon.h" 
#include "Edge.h" 
#include "Vertex.h"
class Cell {
private:
    int id_;
    double volume_;
    double area_;
    double Kv_;
    double V0_;
    double A0_;
    std::array<double, 3> centroid_;
    std::vector<Polygon*> polygons_;
    std::vector<Vertex*> vertices_; 
    void updateVolume();
    void updateArea();
    void updateCentroid();
    void checkPolygonOrientations();
public:
    explicit Cell(std::vector<Vertex*> vertices, std::vector<Polygon*> polygons, int id, double V0, double A0);
    const double getVolume() const;
    const double getArea() const; 
    const double getKv() const;
    const double getV0() const;
    const double getA0() const;
    const std::array<double, 3>& getCentroid() const;
    const int getId() const;
    const std::vector<Polygon*>& getPolygons() const;
    const std::vector<Vertex*>& getVertices() const;
    void update();
    void setId(int id);
    void setKv(double Kv);
    void setV0(double V0);
    void setA0(double A0);
};
