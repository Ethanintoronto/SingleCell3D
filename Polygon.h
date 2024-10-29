#pragma once
#include "Vertex.h"
#include "Edge.h"
#include <vector>
class Polygon{
private:
    int id_;

    double area_;

    double perimeter_;

    std::vector<Vertex*> vertices_;

    std::vector<Edge*> edges_;

    std::array<double, 3> polygonForce_;

    std::array<double, 3> areaVector_;

    std::array<double, 3> centroid_;
    
    void updateCentroid();

    void updateArea();

    void updatePerimeter();

public:
    explicit Polygon(std::vector<Vertex*> vertices, std::vector<Edge*> edges, int id);

    const int getId() const;

    const double getArea() const;

    const double getPerimeter() const;

    const std::vector<Edge*>& getEdges() const;

    const std::vector<Vertex*>& getVertices() const;

    const std::array<double, 3>& getAreaVector() const;

    const std::array<double, 3>& getCentroid() const;

    void setId(int id);

    void update();
    
};