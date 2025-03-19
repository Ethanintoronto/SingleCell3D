#include "Polygon.h"
#include <cmath>
#include <iostream>
Polygon::Polygon(std::vector<Vertex*> vertices, std::vector<Edge*> edges, int id): vertices_(vertices), edges_(edges), id_(id), Ks_(1.), areaVector_({0.,0.,0.}), centroid_({0.,0.,0.}), area_(0.), perimeter_(0.), gamma_(0.),boundary_(false){
    update();
}

int Polygon::getId() const{
    return id_;
}

double Polygon::getArea() const{
    return area_;
}

double Polygon::getPerimeter() const{
    return perimeter_;
}

bool Polygon::isBoundary() const{
    return boundary_;
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

double Polygon::getKs() const{
    return Ks_;
}

double Polygon::getGamma() const{
    return gamma_;
}

void Polygon::setId(int id){
    id_ = id;
}

void Polygon::setKs(double Ks){
    Ks_ = Ks;
}

void Polygon::setGamma(double gamma){
    gamma_ = gamma;
}

void Polygon::setBoundary(bool boundary){
    boundary_ = boundary;
}

void Polygon::updateCentroid(){
    //the following method of calculating the centroid is inspired by ZhangTao's tvm code and uses a method of geometric decomposition
    
    const auto& tmpOrigin = edges_[0]->getVertices()[0]->getPos(); //reference vertex is the first vertex of the first edge in the polygon.
    std::array<double, 3> sumL = {0.,0.,0.};
    double totalLength = 0.;
    for (const auto& edge : edges_){
        double length = edge->getLength();
        const auto& center = edge->getCentroid();
        for (int i = 0; i<3; i++){
            sumL[i]+= length * (center[i]-tmpOrigin[i]);
        }
        totalLength += length; 
    }
    for (int i = 0; i<3; i++){
        centroid_[i] = sumL[i]/totalLength+tmpOrigin[i];
    }

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
    area_ = 0.;
    areaVector_ = {0., 0., 0.};
    std::array<double, 3> ci, cj;
    for (const auto& edge : edges_){
        const auto& pos1 = edge->getVertices()[0]->getPos();
        const auto& pos2 = edge->getVertices()[1]->getPos();

        //Form two vectors from the polygon center to the edge vertices. 
        for (int i = 0; i<3; i++){
            ci[i] = pos1[i]-centroid_[i];
            cj[i] = pos2[i]-centroid_[i];
        }
        //Take the cross product of the two vectors to get the parallelogram area vector  
        std::array<double, 3> crossProduct = {
            ci[1]*cj[2]-ci[2]*cj[1], 
            ci[2]*cj[0]-ci[0]*cj[2],
            ci[0]*cj[1]-ci[1]*cj[0]
        };
        for (int i = 0; i<3; i++){
            areaVector_[i]+=crossProduct[i];
        }
        //Take half the norm of the area vector to get the triangular area
        area_ += 0.5*std::sqrt(crossProduct[0] * crossProduct[0] + 
                               crossProduct[1] * crossProduct[1] + 
                               crossProduct[2] * crossProduct[2]);
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

