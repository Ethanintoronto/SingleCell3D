#include "Cell.h"
#include <cmath>
#include <iostream>
#include <cmath>
Cell::Cell(std::vector<Vertex*> vertices, std::vector<Polygon*> polygons, int id, double V0, double A0): vertices_(vertices), polygons_(polygons), id_(id), V0_(V0), A0_(A0), Kv_(10.),Ka_(1.), area_(0.), centroid_({0.,0.,0.}),geoCentroid_({0.,0.,0.}), volume_(0.){
    update();
    checkPolygonOrientations();
}

int Cell::getId() const{
    return id_;
}

double Cell::getArea() const{
    return area_;
}

double Cell::getVolume() const{
    return volume_;
}

double Cell::getKv() const{
    return Kv_;
}

double Cell::getKa() const{
    return Ka_;
}

double Cell::getV0() const{
    return V0_;
}

double Cell::getA0() const{
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

const std::array<double, 3>& Cell::getGeoCentroid() const{
    return geoCentroid_;
}

double Cell::getEnergy() const{
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

void Cell::updateGeometry(){
    //Updates the centroid, area and volume of the cell 
    
    //the following method is the simplest way of calculating a centroid from the average of the vertex positions:
    double nV = 0; 
    std::array<double, 3> positionSum = {0.,0.,0.};
    for (const auto& vert : vertices_){
        nV++;
        const auto& pos = vert->getPos();
        for (int i =0; i<3; i++){
            positionSum[i]+= pos[i];
        }
    }
    for (int i = 0; i<3; i++){
        centroid_[i] = positionSum[i]/nV;
    }
    
    //The following method computes the volume of the cell and  
    //The weighted geometric centroid from tetrahedron decomposition

    for (int i = 0; i<3; i++){
        geoCentroid_[i] = 0.;
    }
    std::array<double, 3> ci; //face center to curr
    std::array<double, 3> cj; //face center to next
    std::array<double, 3> cc; //cell center to face center
    std::array<double, 3> ctetra; // center of the tetrahedron

    double newVolume = 0.; 
    area_ = 0.;
    for (int i=0; i<polygons_.size(); i++){
        Polygon* polygon = polygons_[i];
        area_ += polygon->getArea();
        for (int j = 0; j<polygon->getVertices().size();j++){
            Vertex* curr = polygon->getVertices()[j];
            Vertex* next = polygon->getVertices()[(j+1)%polygon->getVertices().size()];
            
            //Form two vectors from the polygon center to the edge vertices. 
            const auto& p1 = curr->getPos();
            const auto& p2 = next->getPos(); 
            const auto& polyCentroid = polygon->getCentroid();

            for (int i = 0; i<3; i++){
                ci[i] = p1[i] - polyCentroid[i];
                cj[i] = p2[i] - polyCentroid[i];
                cc[i] = polyCentroid[i] - centroid_[i];
                ctetra[i] = (p1[i]+p2[i]+polyCentroid[i]+centroid_[i])/4;
            }
            
            //Take the cross product of curr cross next to get vector pointing out the cell 
            //and dot with vector from cell center to face center to get positive signed volume
            
            double tetraVolume = ((ci[1]*cj[2]-ci[2]*cj[1])*cc[0] + (ci[2]*cj[0]-ci[0]*cj[2])*cc[1] + (ci[0]*cj[1]-ci[1]*cj[0])*cc[2])/6; 

            geoCentroid_[0] += ctetra[0]*tetraVolume;
            geoCentroid_[1] += ctetra[1]*tetraVolume;
            geoCentroid_[2] += ctetra[2]*tetraVolume;

            newVolume += tetraVolume;
        }
    }
    volume_= newVolume;

    for (int i =0; i<3;i++){
        geoCentroid_[i] /= volume_;
    }     
}

void Cell::update(){
    for (auto& polygon : polygons_){
        polygon->update();    
    }
    updateGeometry();
}

void Cell::checkPolygonOrientations() {
        for (const auto& polygon : polygons_) {
            const auto& pv = polygon->getVertices();
            const auto& polyCentroid = polygon->getCentroid(); 
            double totalCrossZ = 0.0;
            std::size_t n = pv.size();
            for (std::size_t i = 0; i < n; ++i) {
                const auto& p1 = pv[i]->getPos();
                const auto& p2 = pv[(i + 1) % n]->getPos();
                std::array<double, 3> ci, cj, cc; 

                for (int i = 0; i<3; i++){
                    ci[i] = p1[i] - polyCentroid[i];
                    cj[i] = p2[i] - polyCentroid[i];
                    cc[i] = polyCentroid[i] - centroid_[i];
                }

                // Cross product z-component
                double crossX = (ci[1] * cj[2]) - (ci[2] * cj[1]);
                double crossY = (ci[2] * cj[0]) - (ci[0] * cj[2]);
                double crossZ = (ci[0] * cj[1]) - (ci[1] * cj[0]);

                // Find the angle between the vector from the cell center to face center and the surface vector
                // a dot b = abcostheta
                // theta = acos(a dot b/|a|b|)
                double cosTheta = (cc[0]*crossX + cc[1]*crossY + cc[2]*crossZ)/(std::sqrt(std::pow(cc[0],2)+std::pow(cc[1],2)+std::pow(cc[2],2))*std::sqrt(std::pow(crossX,2)+std::pow(crossY,2)+std::pow(crossZ,2)));
                //clamp values to [-1, 1]
                cosTheta = std::max(-1.0, std::min(1.0, cosTheta));
                double theta = std::acos(cosTheta);
                if (theta>=M_PI_2){
                    throw std::invalid_argument("Polygon orientation error");
                }
            }
        }
    }
