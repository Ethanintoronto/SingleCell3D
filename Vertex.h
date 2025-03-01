#pragma once 
#include <array>
#include <vector>
class Vertex{
private:
    std::array<double,3> pos_; //Array of 3D position coordinates 
    std::array<double,3> force_; //Force acting on the vertex
    std::vector<std::array<double,3>> posHistory_; 
    std::vector<std::array<double, 3>> forceHistory_;
    int id_; //Note if the number of vertices in a simulation exceeds 2,147,483,647 long int datatype should be used
public:
    explicit Vertex(std::array<double,3> pos, int id); //Constructor for vertex with initial position
    
    const std::array<double,3>& getPos() const; //Getter method for the vertex position 

    const std::array<double,3>& getForce() const; //Getter method for the vertex position 
    
    int getId() const; //Getter method for the vertex id

    void setPos(std::array<double,3> pos); //Setter method for the vertex position

    void setForce(std::array<double, 3> force); //Setter method for the vertex force

    void setId(int id);

    void updateHist();
};