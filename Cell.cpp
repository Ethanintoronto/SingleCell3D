#include "Cell.h"
#include <cmath>
Cell::Cell(std::vector<Vertex*> vertices, std::vector<Polygon*> polygons, int id): vertices_(vertices), polygons_(polygons), id_(id), area_(0.), centroid_({0.,0.,0.}), volume_(0.){
    update();
}

const int Cell::getId() const{
    return id_;
}

const double Cell::getArea() const{
    return area_;
}

const double Cell::getArea() const{
    return area_;
}

const double Cell::getVolume() const{
    return volume_;
}
const std::vector<Vertex*>& Cell::getVertices() const{
    return vertices_;
}

const std::vector<Polygon*>& Cell::getPolygons() const{
    return polygons_;
}

const std::array<double, 3>& Cell::getCentroid() const{
    return centroid_;
}

void Cell::setId(int id){
    id_ = id;
}

void Cell::updateCentroid(){
    //the following method is the simplest way of calculating a centroid from the average vertex positions:
    double nV = 0; 
    double dx[3];
    for (const auto& vert : vertices_){
        nV++;
        std::array<double, 3> pos = vert->getPos();
        dx[0] += pos[0];
        dx[1] += pos[1];
        dx[2] += pos[2];
    }
    centroid_[0]=dx[0]/nV;
    centroid_[1]=dx[1]/nV;
    centroid_[2]=dx[2]/nV;
}

void Cell::updateVolume(){
    //the following method uses a triangulation scheme taking the tetrahedron volumes formed by the cell centroid two edge vertices and the face centroids.
    double newVolume = 0; 
    for (const auto& polygon : polygons_){
        std::array<double,3> polyCenter = polygon->getCentroid(); 
        double ci[3] = {0,0,0};
        double cj[3] = {0,0,0};
        double cc[3] = {0,0,0};
        for (const auto& edge : polygon->getEdges()){
            std::array<double, 3> pos1 = edge->getVertices()[0]->getPos();
            std::array<double, 3> pos2 = edge->getVertices()[1]->getPos();
            ci[0] = pos1[0]-centroid_[0];
            cj[0] = pos2[0]-centroid_[0];
            ci[1] = pos1[1]-centroid_[1];
            cj[1] = pos2[1]-centroid_[1];
            ci[2] = pos1[2]-centroid_[2];
            cj[2] = pos2[2]-centroid_[2];
            cc[0] = polyCenter[0]-centroid_[0];
            cc[1] = polyCenter[1]-centroid_[1];
            cc[2] = polyCenter[2]-centroid_[2];
        newVolume += std::abs((ci[1]*cj[2]-ci[2]*cj[1])*cc[0]);
        newVolume += std::abs((ci[2]*cj[0]-ci[0]*cj[2])*cc[1]);
        newVolume += std::abs((ci[0]*cj[1]-ci[1]*cj[0])*cc[2]);
        }
    }
    volume_ = newVolume/6;
}
void Cell::updateArea(){
    double newArea = 0; 
    for (const auto& polygon : polygons_){
        newArea+=polygon->getArea();
    area_ = newArea;
    }
}

void Cell::update(){
    for (auto& polygon : polygons_){
        polygon->update();    
    }
    updateCentroid();
    updateArea();
    updateVolume();
}


