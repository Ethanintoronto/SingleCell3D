#pragma once
#include "Cell.h"
class Simulation{
private:
    double time_;
    double timestep_;
    int numTimesteps_; 
    double Kv_;
    double Ka_;
    double V0_;
    double A0_;
    double eta_;
    int log_;
    std::vector<Cell*> cells_; 
    std::vector<Vertex*> vertices_;
    void Run();
    void updateForces();
    void updateCells();
    void performTimeStep();
    void update();
    std::array<double, 3> dAdr(Vertex* current, Vertex* prev, Vertex* next, std::array<double,3> polyCenter, std::array<double,3> cellCenter, int N_p);
    std::array<double, 3> dVdr(Vertex* current, Vertex* prev, Vertex* next, std::array<double,3> polyCenter, std::array<double,3> cellCenter, int N_p, int N_c);
public:
    explicit Simulation(std::vector<Cell*> cells, std::vector<Vertex*> vertices, double timestep, int numTimesteps, double Kv, double Ka, double V0, double A0, double eta,int log);
};