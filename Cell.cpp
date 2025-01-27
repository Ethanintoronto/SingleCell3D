#include "Cell.h"
#include <cmath>
#include <iostream>
Cell::Cell(std::vector<Vertex*> vertices, std::vector<Polygon*> polygons, int id, double V0, double A0): vertices_(vertices), polygons_(polygons), id_(id), V0_(V0), A0_(A0), Kv_(10.),Ka_(1.), area_(0.), centroid_({0.,0.,0.}), volume_(0.){
    update();
    checkPolygonOrientations();
}

const int Cell::getId() const{
    return id_;
}

const double Cell::getArea() const{
    return area_;
}

const double Cell::getVolume() const{
    return volume_;
}

const double Cell::getKv() const{
    return Kv_;
}

const double Cell::getKa() const{
    return Ka_;
}

const double Cell::getV0() const{
    return V0_;
}

const double Cell::getA0() const{
    return A0_;
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

const double Cell::getEnergy() const{
    return Ka_*std::pow((area_-A0_),2) + Kv_*std::pow((volume_-V0_),2);
}

void Cell::setId(int id){
    id_ = id;
}

void Cell::setKv(double Kv){
    Kv_ = Kv;
}

void Cell::setKa(double Ka){
    Ka_ = Ka;
}

void Cell::setV0(double V0){
    V0_ = V0;
}

void Cell::setA0(double A0){
    A0_ = A0;
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
        double ci[3] = {0,0,0}; //vector from the cell centroid to the first edge vertex 
        double cj[3] = {0,0,0}; //vector from the cell centroid to the second edge vertex
        double cc[3] = {0,0,0}; //vector from the cell centroid to the polygon center
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
        //Taking the triple product of the three vectors to compute the volume of the
        //tetrahedron with vertices at the cell center, face center, and edges.
        newVolume += std::abs((ci[1]*cj[2]-ci[2]*cj[1])*cc[0]);
        newVolume += std::abs((ci[2]*cj[0]-ci[0]*cj[2])*cc[1]);
        newVolume += std::abs((ci[0]*cj[1]-ci[1]*cj[0])*cc[2]);
        }
    }
    volume_ = newVolume/6;
}
void Cell::updateArea(){
    area_ = 0.; 
    for (const auto& polygon : polygons_){
        area_+=polygon->getArea();
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

void Cell::checkPolygonOrientations() {
        for (const auto& polygon : polygons_) {
            const auto& pv = polygon->getVertices();
            double totalCrossZ = 0.0;
            std::size_t n = pv.size();
            for (std::size_t i = 0; i < n; ++i) {
                const auto& v1 = pv[i]->getPos();
                const auto& v2 = pv[(i + 1) % n]->getPos();

                // Vectors from centroid to vertices
                std::array<double, 3> vec1 = {v1[0] - polygon->getCentroid()[0], v1[1] - polygon->getCentroid()[1], v1[2] - polygon->getCentroid()[2]};
                std::array<double, 3> vec2 = {v2[0] - polygon->getCentroid()[0], v2[1] - polygon->getCentroid()[1], v2[2] - polygon->getCentroid()[2]};

                // Cross product z-component
                double crossX = (vec1[1] * vec2[2]) - (vec1[2] * vec2[1]);
                double crossY = (vec1[2] * vec2[0]) - (vec1[0] * vec2[2]);
                double crossZ = (vec1[0] * vec2[1]) - (vec1[1] * vec2[0]);

                // Vector from cell center to face center:
                std::array<double, 3> centerFace = {polygon->getCentroid()[0]-centroid_[0], polygon->getCentroid()[1]-centroid_[1], polygon->getCentroid()[2]-centroid_[2]};
                // Find the angle between the vector from the cell center to face center and the surface vector
                // a dot b = abcostheta
                // theta = acos(a dot b/|a|b|)
                double cosTheta = (centerFace[0]*crossX + centerFace[1]*crossY + centerFace[2]*crossZ)/(std::sqrt(std::pow(centerFace[0],2)+std::pow(centerFace[1],2)+std::pow(centerFace[2],2))*std::sqrt(std::pow(crossX,2)+std::pow(crossY,2)+std::pow(crossZ,2)));
                //clamp values to [-1, 1]
                cosTheta = std::max(-1.0, std::min(1.0, cosTheta));
                double theta = std::acos(cosTheta);
                if (theta>=1.57079632679){
                    std::cout<<"WARNING: Polygon orientation error";
                }
            }
        }
    }


