#include "Polygon.h"
#include <cmath>
#include <iostream>
Polygon::Polygon(std::vector<Vertex*> vertices, std::vector<Edge*> edges, int id): vertices_(vertices), edges_(edges), id_(id), polygonForce_({0.,0.,0.}), areaVector_({0.,0.,0.}), centroid_({0.,0.,0.}), area_(0.), perimeter_(0.){
    update();
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

const std::vector<Vertex*>& Polygon::getVertices() const{
    return vertices_;
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

    //the following method is the simplest way of calculating a centroid from the average vertex positions:
    /*
    double dx[3];
    for (const auto& vert : vertices_){
        std::array<double, 3> pos = vert->getPos();
        dx[0] += pos[0];
        dx[1] += pos[1];
        dx[2] += pos[2];
    }
    centroid_[0]=dx[0]/vertices_.size();
    centroid_[1]=dx[1]/vertices_.size();
    centroid_[2]=dx[2]/vertices_.size();
    */
}

void Polygon::updateArea(){
    double areaVectorX;
    double areaVectorY;
    double areaVectorZ;
    area_ = 0.;
    areaVector_[0] = 0.;
    areaVector_[1] = 0.;
    areaVector_[2] = 0.;
    std::array<double, 3> ci;
    std::array<double, 3> cj;
    for (const auto& edge : edges_){
        std::array<double, 3> pos1 = edge->getVertices()[0]->getPos();
        std::array<double, 3> pos2 = edge->getVertices()[1]->getPos();
        ci[0] = pos1[0]-centroid_[0];
        cj[0] = pos2[0]-centroid_[0];

        ci[1] = pos1[1]-centroid_[1];
        cj[1] = pos2[1]-centroid_[1];

        ci[2] = pos1[2]-centroid_[2];
        cj[2] = pos2[2]-centroid_[2];
        areaVectorX = (ci[1]*cj[2]-ci[2]*cj[1]);
        areaVectorY = (ci[2]*cj[0]-ci[0]*cj[2]);
        areaVectorZ = (ci[0]*cj[1]-ci[1]*cj[0]);

        areaVector_[0] += areaVectorX;
        areaVector_[1] += areaVectorY;
        areaVector_[2] += areaVectorZ;
        area_ += 0.5*std::sqrt(std::pow(areaVectorX,2)+std::pow(areaVectorY,2)+std::pow(areaVectorZ,2));
    }
}

void Polygon::updatePerimeter(){
    double newPerimeter = 0.;
    for (const auto& edge : edges_){
        newPerimeter+= edge->getLength();
    }
    perimeter_ = newPerimeter;
}

void Polygon::update(){
    for (auto& edge : edges_){
        edge->update();    
    }

    updateCentroid();
    updateArea();
    updatePerimeter();
}

