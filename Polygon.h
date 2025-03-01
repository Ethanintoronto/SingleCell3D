#pragma once
#include "Vertex.h"
#include "Edge.h"
class Polygon{
private:
    int id_;

    double area_;

    double perimeter_;

    double Ka_;

    double gamma_;

    std::vector<Vertex*> vertices_;

    std::vector<Edge*> edges_;

    std::array<double, 3> areaVector_;

    std::array<double, 3> centroid_;
    
    void updateCentroid();

    void updateArea();

    void updatePerimeter();

public:
    explicit Polygon(std::vector<Vertex*> vertices, std::vector<Edge*> edges, int id);

    int getId() const;

    double getArea() const;

    double getPerimeter() const;

    double getKa() const;

    const std::vector<Edge*>& getEdges() const;

    const std::vector<Vertex*>& getVertices() const;

    const std::array<double, 3>& getAreaVector() const;

    const std::array<double, 3>& getCentroid() const;

    double getGamma() const;

    void setGamma(double gamma);

    void setId(int id);

    void setKa(double Ka);

    void update();
    
};