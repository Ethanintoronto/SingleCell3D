#include "Edge.h"
#include <cmath>

Edge::Edge(std::array<Vertex*,2> vertices, int id): vertices_(vertices), id_(id) {
    update();
}

const std::array<Vertex*,2>& Edge::getVertices() const{
    return vertices_;
}

const double Edge::getLength() const{
    return length_;
}

const int Edge::getId() const{
    return id_;
}

const std::array<double, 3>& Edge::getCentroid() const{
    return centroid_;
}

void Edge::setId(int id){
    id_ = id;
}

void Edge::updateLength() {
    const auto& pos1 = vertices_[0]->getPos(); // Using const reference to avoid copying
    const auto& pos2 = vertices_[1]->getPos(); // Using const reference to avoid copying

    // Calculate the squared differences
    double dx = pos2[0] - pos1[0];
    double dy = pos2[1] - pos1[1];
    double dz = pos2[2] - pos1[2];

    // Calculate length using the squared differences
    length_ = std::sqrt(dx * dx + dy * dy + dz * dz);
}

void Edge::updateCentroid(){
    const auto& pos1 = vertices_[0]->getPos();
    const auto& pos2 = vertices_[1]->getPos();

    //Calculate centroid using the average position of the two vertices:
    //Must change to implement periodic boundary conditions
    for (int i = 0; i<3;i++){
        centroid_[i] = (pos2[i]+pos1[i])/2;
    }
}

void Edge::update(){
    updateCentroid();
    updateLength();
}


