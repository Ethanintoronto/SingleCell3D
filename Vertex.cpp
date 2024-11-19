#include "Vertex.h"
Vertex::Vertex(std::array<double, 3> pos, int id): pos_(pos), force_({0.,0.,0.}), id_(id){}

const std::array<double,3>& Vertex::getPos() const{
    return pos_;
}

const std::array<double,3>& Vertex::getForce() const{
    return force_;
}

const int Vertex::getId() const{
    return id_;
}

void Vertex::setPos(std::array<double, 3> pos){
    pos_ = pos;
}
void Vertex::setForce(std::array<double, 3> force){
    force_ = force;
}

void Vertex::setId(int id){
    id_ = id;
}
void Vertex::updateHist(){
    posHistory_.push_back(pos_);
    forceHistory_.push_back(force_);
}


