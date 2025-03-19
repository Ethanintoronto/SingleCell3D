#include "Cell.h"
#include <cmath>
#include <iostream>
Cell::Cell(std::vector<Vertex*> vertices, std::vector<Polygon*> polygons, int id, double V0, double A0): vertices_(vertices), polygons_(polygons), id_(id), V0_(V0), A0_(A0), Kv_(10.),Ka_(1.), area_(0.), centroid_({0.,0.,0.}),geoCentroid_({0.,0.,0.}), volume_(0.){
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
    checkPolygonOrientations();
    update();
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

const std::vector<bool>& Cell::getPolyOrientations() const{
    return polyOrientations_;
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
    
    //The following method computes the volume of the cell and  
    //The weighted geometric centroid from tetrahedron decomposition

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
            Vertex* curr;
            Vertex* next;
            if (polyOrientations_[i]){
                curr = polygon->getVertices()[j];
                next = polygon->getVertices()[(j+1)%polygon->getVertices().size()];
            }
            else{
                curr = polygon->getVertices()[j];
                if (j==0){
                    next = polygon->getVertices()[polygon->getVertices().size()-1];
                }
                else{
                    next = polygon->getVertices()[j-1];
                }
            }
            
            //Form two vectors from the polygon center to the edge vertices. 
            const auto& p1 = curr->getPos();
            const auto& p2 = next->getPos(); 
            const auto& polyCentroid = polygon->getCentroid();

            for (int k = 0; k<3; k++){
                ci[k] = p1[k] - polyCentroid[k];
                cj[k] = p2[k] - polyCentroid[k];
                cc[k] = polyCentroid[k] - centroid_[k];
                ctetra[k] = (p1[k]+p2[k]+polyCentroid[k]+centroid_[k])/4;
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
        for (int h =0;h<polygons_.size(); h++) {
            const auto& polygon = polygons_[h]; 
            const auto& pv = polygon->getVertices();
            const auto& polyCentroid = polygon->getCentroid(); 
            std::array<double, 3> cc;
            for (int j=0;j<3;j++){
                cc[j] = polyCentroid[j] - centroid_[j];
            }
            
            for (int i = 0; i < pv.size(); i++) {
                const auto& p1 = pv[i]->getPos();
                const auto& p2 = pv[(i + 1) % pv.size()]->getPos();
                std::array<double, 3> ci, cj;

                for (int j = 0; j<3; j++){
                    ci[j] = p1[j] - polyCentroid[j];
                    cj[j] = p2[j] - polyCentroid[j];
                }
               

                // Cross product z-component
                double crossX = (ci[1] * cj[2]) - (ci[2] * cj[1]);
                double crossY = (ci[2] * cj[0]) - (ci[0] * cj[2]);
                double crossZ = (ci[0] * cj[1]) - (ci[1] * cj[0]);

                // Find the angle between the vector from the cell center to face center and the surface vector
                // a dot b = |a| |b| costheta
                // theta = acos(a dot b/|a|b|)
                double cosTheta = (cc[0]*crossX + cc[1]*crossY + cc[2]*crossZ)/(std::sqrt(cc[0]*cc[0]+cc[1]*cc[1]+cc[2]*cc[2])*std::sqrt(crossX*crossX+crossY*crossY+crossZ*crossZ));
                //clamp values to [-1, 1]
                cosTheta = std::max(-1.0, std::min(1.0, cosTheta));

                double theta = std::acos(cosTheta);
                
                if (i==0){
                    if (theta>=M_PI_2){
                        //throw std::invalid_argument("Polygon orientation error");
                        polyOrientations_.push_back(false);
                    }
                    else {
                        polyOrientations_.push_back(true);
                    }   
                }
                else if ((theta>=M_PI_2 && polyOrientations_[h])||(theta<M_PI_2&&!polyOrientations_[h])){
                    throw std::runtime_error("Inconsistent polygon direction in Cell: "+std::to_string(id_)+" Polygon: "+std::to_string(polygon->getId()) + " Vertex: " + std::to_string(i));
                }
            }
        }
    }
