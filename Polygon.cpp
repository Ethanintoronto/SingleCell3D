#include "Polygon.h"

Polygon::Polygon(std::vector<Edge*> edges, int id): edges_(edges), id_(id), polygonForce_({0.,0.,0.}), areaVector_({0.,0.,0.}), centroid_({0.,0.,0.}){
    updateArea();
}

const int Polygon::getId() const{
    return id_;
}

const double Polygon::getArea() const{
    return area_;
}

const double Polygon::getPerimeter() const{
    return perimeter_;
}

const std::vector<Edge*>& Polygon::getEdges() const{
    return edges_;
}

const std::array<double, 3>& Polygon::getAreaVector() const{
    return areaVector_;
}

const std::array<double, 3>& Polygon::getCentroid() const{
    return centroid_;
}

void Polygon::setId(int id){
    id_ = id;
}

void Polygon::updateCentroid(){
    //the following method of calculating the centroid is inspired by ZhangTao's tvm code and uses a method of geometric decomposition
    const std::array<double, 3>& tmpOrigin = edges_[0]->getVertices()[0]->getPos(); //reference vertex is the first vertex of the first edge in the polygon.
    double sumLx = 0.;
    double sumLy = 0.;
    double sumLz = 0.;
    double sumL = 0.;
    for (const auto& edge : edges_){
        double length = edge->getLength();
        double dx[3];
        const auto& center = edge->getCentroid();
        dx[0] = center[0] - tmpOrigin[0];
        dx[1] = center[1] - tmpOrigin[1];
        dx[2] = center[2] - tmpOrigin[2];
        //check periodic BC on dx here
        sumLx+=length*dx[0];
        sumLy+=length*dx[1];
        sumLz+=length*dx[2];
        sumL += length; 
    }
    centroid_[0] = sumLx/sumL + tmpOrigin[0];
    centroid_[1] = sumLy/sumL + tmpOrigin[1];
    centroid_[2] = sumLz/sumL + tmpOrigin[2];
}

void Polygon::updateArea(){

}

void Polygon::update(){
    for (auto& edge : edges_){
        edge->update();    
    }

    updateCentroid();
    updateArea();
}

