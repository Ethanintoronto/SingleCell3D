#pragma once
#include "Cell.h"
#include <string>
class Simulation{
private:
    int time_;
    int id_;
    double timestep_;
    int numTimesteps_;
    double eta_;
    int log_;
    bool write_;
    double period_;
    std::vector<Cell*> cells_; 
    std::vector<Vertex*> vertices_;
    std::vector<Edge*> edges_;
    std::vector<Polygon*> polygons_;
    void Run();
    void updateForces();
    void updateCells();
    void performTimeStep();
    void writeVTK();
    void printMaxForce();
    void writeMaxForce();
    void writeVolume();
    void writeArea();
    void writeEnergy(); 
    void writeCellCentroid();
    void writeCellGeoCentroid();
    void update();
    void setId(int id);
    void setPeriod(double period);
    std::string convertDouble(double val);
    std::string getDate();
    const int getId() const;
    const double getPeriod() const;
    std::array<double, 3> dAdr(Vertex* current, Vertex* prev, Vertex* next, std::array<double,3> polyCenter, int N_p);
    std::array<double, 3> dVdr(Vertex* current, Vertex* prev, Vertex* next, std::array<double,3> polyCenter, std::array<double,3> cellCenter, int N_p, int N_c);
    std::array<double, 3> dAShoeLace(Vertex* prev, Vertex* next, int faceNorm);

public:
    explicit Simulation(std::vector<Cell*> cells, std::vector<Polygon*> polygons, std::vector<Edge*> edges, std::vector<Vertex*> vertices, int id, double period, double timestep, int numTimesteps, double eta, int log, bool write);
};