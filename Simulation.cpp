#include "Simulation.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <algorithm>
Simulation::Simulation(std::vector<Cell*> cells, std::vector<Polygon*> polygons, std::vector<Edge*> edges, std::vector<Vertex*> vertices, double timestep, int numTimesteps, double Kv, double Ka, double V0, double A0, double eta, int log, bool write): cells_(cells), polygons_(polygons), edges_(edges), vertices_(vertices), timestep_(timestep),numTimesteps_(numTimesteps), Kv_(Kv), Ka_(Ka), V0_(V0), A0_(A0), eta_(eta), log_(log),time_(0),write_(write){
    Run();
}
void Simulation::Run(){
    if (write_ && time_%log_==0){
        writeVTK();
        writeMaxForce();
    }

    for (int i=0;i<numTimesteps_;i++){
        //std::cout<<"Cell Area: "<<cells_[0]->getArea()<<"\n";
        //std::cout<<"Cell Volume: "<<cells_[0]->getVolume()<<"\n";
        update();
        if (write_ && time_%log_==0){
            writeVolume();
            writeArea();
            writeCellCentroid();
        }
        //printMaxForce();
        //std::cout<<"Cell Centroid: ("<<cells_[0]->getCentroid()[0]<<","<<cells_[0]->getCentroid()[1]<<","<<cells_[0]->getCentroid()[2]<<")\n";
        time_ += 1;
        if (write_ && time_%log_==0){
            writeVTK();
            writeMaxForce();
        }
    }
}
void Simulation::updateForces(){
    //Reset all forces
    for (auto& v :vertices_){
        v->setForce({0.,0.,0.});
    }
    int prev_k;
    int next_k; 
    std::array<double, 3> pos;
    for (int i = 0; i<cells_.size();i++){
        for (int j = 0;  j<cells_[i]->getPolygons().size(); j++){
            Polygon* polygon = cells_[i]->getPolygons()[j];
            for (int k =  0; k<polygon->getVertices().size();k++){
                Vertex* curr = polygon->getVertices()[k];
                std::array<double,3> force =  curr->getForce();
                if (k==0){
                    prev_k = polygon->getVertices().size()-1;
                    next_k = k+1;
                }
                else if (k==polygon->getVertices().size()-1){
                    prev_k = k-1;
                    next_k = 0;
                }
                else{
                    prev_k = k-1;
                    next_k = k+1; 
                }
                std::array<double,3> areaDeriv = dAdr(curr,polygon->getVertices()[prev_k],polygon->getVertices()[next_k],polygon->getCentroid(),cells_[i]->getCentroid(),polygon->getVertices().size());
                std::array<double,3> volDeriv = dVdr(curr, polygon->getVertices()[prev_k], polygon->getVertices()[next_k],polygon->getCentroid(), cells_[i]->getCentroid(), polygon->getVertices().size(), cells_[i]->getVertices().size());
                for (int l =0; l<3; l++){
                    force[l] -= 2*Ka_*(cells_[i]->getArea()-A0_)*areaDeriv[l] + 2*Kv_*(cells_[i]->getVolume()-V0_)*volDeriv[l];
                }
                polygon->getVertices()[k]->setForce(force);
            }
        }
    }
}

void Simulation::updateCells(){
    for (const auto& cell : cells_){
        cell->update();
    }
}

void Simulation::update(){
    updateCells();
    updateForces();
    performTimeStep();
}

void Simulation::performTimeStep(){
    for (const auto& vertex : vertices_){
        vertex->updateHist();
        std::array<double, 3> newpos = {0,0,0};
        for (int i=0; i<3; i++){
            newpos[i] = vertex->getPos()[i] + eta_*vertex->getForce()[i]*timestep_; 
        }
        vertex->setPos(newpos);
    }
}

// Function to write simulation data in VTK format
void Simulation::writeVTK() {
    std::ostringstream dataDir;
    dataDir << "data/vtk/Single_cell_sim_V0_"<<convertDouble(V0_)<<"_A0_"<<convertDouble(A0_)<<"_timestep_"<<convertDouble(timestep_);
    std::cout << dataDir.str()<<"\n";
    if (!std::filesystem::exists(dataDir.str())){
        if(std::filesystem::create_directory(dataDir.str())){
            std::cout <<"Datadir created successfully\n";
        }
        else{
            std::cout <<"Failed to create directory\n";
        }
    }
    std::ostringstream filename;
    filename << dataDir.str()<<"/Single_cell_sim_V0_"<<convertDouble(V0_)<<"_A0_"<<convertDouble(A0_)<<"_timestep_"<<convertDouble(timestep_)<<"_"<< std::setw(3) << std::setfill('0') << time_ << ".vtk";
    std::ofstream vtkFile(filename.str());

    if (!vtkFile.is_open()) {
        std::cerr << "Error opening the VTK file for writing!" << std::endl;
        return;
    }

    // Write the header for VTK file
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Simulation Data at timestep " << time_ << "\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET POLYDATA\n";

    // Write vertices
    vtkFile << "POINTS " << vertices_.size() << " float\n";
    for (const auto& vertex : vertices_) {
        const auto& pos = vertex->getPos();
        vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
    }

    // Write edges
    vtkFile << "LINES " << edges_.size() << " " << edges_.size() * 3 << "\n";
    for (const auto& edge : edges_) {
        vtkFile << "2 " << edge->getVertices()[0]->getId() << " "
                << edge->getVertices()[1]->getId() << "\n";
    }

    // Write polygons
    vtkFile << "POLYGONS " << polygons_.size() << " " << polygons_.size() * 5 << "\n";
    for (const auto& polygon : polygons_) {
        vtkFile << polygon->getVertices().size();
        for (const auto& vertex : polygon->getVertices()) {
            vtkFile << " " << vertex->getId();
        }
        vtkFile << "\n";
    }

    // Close the VTK file
    vtkFile.close();
    std::cout << "VTK file for timestep " << time_ << " has been written successfully.\n";
}

void Simulation::printMaxForce(){
    // Find the max force:
    double max_force = 0.;
    for (auto& vertex : vertices_){
        if (std::sqrt(std::pow(vertex->getForce()[0],2)+std::pow(vertex->getForce()[1],2)+std::pow(vertex->getForce()[2],2)) > max_force){
            max_force = std::sqrt(std::pow(vertex->getForce()[0],2)+std::pow(vertex->getForce()[1],2)+std::pow(vertex->getForce()[2],2)); 
        }
    }
    std::cout<<"Max force: "<<max_force<<"\n";
}

void Simulation::writeMaxForce(){
    // Find the max force:
    double max_force = 0.;
    int max_id = 0;
    int curr_id = 0;
    for (auto& vertex : vertices_){
        if (std::sqrt(std::pow(vertex->getForce()[0],2)+std::pow(vertex->getForce()[1],2)+std::pow(vertex->getForce()[2],2)) > max_force){
            max_force = std::sqrt(std::pow(vertex->getForce()[0],2)+std::pow(vertex->getForce()[1],2)+std::pow(vertex->getForce()[2],2)); 
            max_id = curr_id;
        }
        curr_id++;
    }

    //Open file:
    std::ostringstream dataDir;
    dataDir << "data/Single_cell_sim_V0_"<<convertDouble(V0_)<<"_A0_"<<convertDouble(A0_)<<"_timestep_"<<convertDouble(timestep_);
    if (!std::filesystem::exists(dataDir.str())){
        if(std::filesystem::create_directory(dataDir.str())){
            std::cout <<"Datadir created successfully\n";
        }
        else{
            std::cout <<"Failed to create directory\n";
        }
    }
    std::ostringstream filename;
    filename << dataDir.str()<<"/Single_cell_MaxForce_V0_"<<convertDouble(V0_)<<"_A0_"<<convertDouble(A0_)<<"_timestep_"<<convertDouble(timestep_)<< ".txt";
    std::ofstream forceFile;
    if (time_ ==0){
        forceFile.open(filename.str(), std::ios::trunc);
    }
    else {
        forceFile.open(filename.str(), std::ios::app);
    }
    if (!forceFile.is_open()) {
        std::cerr << "Error opening the force file for writing!" << std::endl;
        return;
    }

    forceFile<<time_*timestep_<<","<< max_force<<","<< vertices_[max_id]->getForce()[0]<<","<< vertices_[max_id]->getForce()[1]<<","<< vertices_[max_id]->getForce()[2]<<"\n";
    // Close the file
    forceFile.close();
}
void Simulation::writeVolume(){
    //Open file:
    std::ostringstream dataDir;
    dataDir << "data/Single_cell_sim_V0_"<<convertDouble(V0_)<<"_A0_"<<convertDouble(A0_)<<"_timestep_"<<convertDouble(timestep_);
    if (!std::filesystem::exists(dataDir.str())){
        if(std::filesystem::create_directory(dataDir.str())){
            std::cout <<"Datadir created successfully\n";
        }
        else{
            std::cout <<"Failed to create directory\n";
        }
    }
    std::ostringstream filename;
    filename << dataDir.str()<<"/Single_cell_Volume_V0_"<<convertDouble(V0_)<<"_A0_"<<convertDouble(A0_)<<"_timestep_"<<convertDouble(timestep_)<< ".txt";
    std::ofstream volumeFile;
    if (time_ ==0){
        volumeFile.open(filename.str(), std::ios::trunc);
    }
    else {
        volumeFile.open(filename.str(), std::ios::app);
    }
    if (!volumeFile.is_open()) {
        std::cerr << "Error opening the force file for writing!" << std::endl;
        return;
    }
    volumeFile<<time_*timestep_<<","<< cells_[0]->getVolume()<<"\n";
    // Close the file
    volumeFile.close();
}

void Simulation::writeArea(){
    //Open file:
    std::ostringstream dataDir;
    dataDir << "data/Single_cell_sim_V0_"<<convertDouble(V0_)<<"_A0_"<<convertDouble(A0_)<<"_timestep_"<<convertDouble(timestep_);
    if (!std::filesystem::exists(dataDir.str())){
        if(std::filesystem::create_directory(dataDir.str())){
            std::cout <<"Datadir created successfully\n";
        }
        else{
            std::cout <<"Failed to create directory\n";
        }
    }
    std::ostringstream filename;
    filename << dataDir.str()<<"/Single_cell_Area_V0_"<<convertDouble(V0_)<<"_A0_"<<convertDouble(A0_)<<"_timestep_"<<convertDouble(timestep_)<< ".txt";
    std::ofstream areaFile;
    if (time_ ==0){
        areaFile.open(filename.str(), std::ios::trunc);
    }
    else {
        areaFile.open(filename.str(), std::ios::app);
    }

    if (!areaFile.is_open()) {
        std::cerr << "Error opening the force file for writing!" << std::endl;
        return;
    }
    areaFile<<time_*timestep_<<","<< cells_[0]->getArea()<<"\n";
    // Close the file
    areaFile.close();
}

void Simulation::writeCellCentroid(){
    //Open file:
    std::ostringstream dataDir;
    dataDir << "data/Single_cell_sim_V0_"<<convertDouble(V0_)<<"_A0_"<<convertDouble(A0_)<<"_timestep_"<<convertDouble(timestep_);
    if (!std::filesystem::exists(dataDir.str())){
        if(std::filesystem::create_directory(dataDir.str())){
            std::cout <<"Datadir created successfully\n";
        }
        else{
            std::cout <<"Failed to create directory\n";
        }
    }
    std::ostringstream filename;
    filename << dataDir.str()<<"/Single_cell_Centroid_V0_"<<convertDouble(V0_)<<"_A0_"<<convertDouble(A0_)<<"_timestep_"<<convertDouble(timestep_)<< ".txt";

    std::ofstream centroidFile;
    if (time_ ==0){
        centroidFile.open(filename.str(), std::ios::trunc);
    }
    else {
        centroidFile.open(filename.str(), std::ios::app);
    }

    if (!centroidFile.is_open()) {
        std::cerr << "Error opening the force file for writing!" << std::endl;
        return;
    }
    centroidFile<<time_*timestep_<<","<< cells_[0]->getCentroid()[0]<<","<< cells_[0]->getCentroid()[1]<<","<< cells_[0]->getCentroid()[2]<<"\n";
    // Close the file
    centroidFile.close();
}


std::array<double, 3> Simulation::dAdr(Vertex* current, Vertex* prev, Vertex* next, std::array<double,3> polyCenter, std::array<double,3> cellCenter, int N_p){
    double x_k = current->getPos()[0];
    double y_k = current->getPos()[1];
    double z_k = current->getPos()[2];

    double x_kp1 = next->getPos()[0];
    double y_kp1 = next->getPos()[1];
    double z_kp1 = next->getPos()[2];

    double x_km1 = prev->getPos()[0];
    double y_km1 = prev->getPos()[1];
    double z_km1 = prev->getPos()[2];

    double Px_0 = polyCenter[0];
    double Py_0 = polyCenter[1];
    double Pz_0 = polyCenter[2];
    double dAdx = (2*x_k*(std::pow(Py_0, 2) - 2*Py_0*y_km1 + std::pow(Pz_0, 2) - 2*Pz_0*z_km1 + std::pow(y_km1, 2) + std::pow(z_km1, 2)) + x_km1*(2*Py_0*y_k + 2*Py_0*y_km1 + 2*Pz_0*z_k + 2*Pz_0*z_km1 - 2*y_k*y_km1 - 2*z_k*z_km1 - 2*std::pow(Py_0, 2) - 2*std::pow(Pz_0, 2)) - Px_0*(2*Py_0*y_k - 2*Py_0*y_km1 + 2*Pz_0*z_k - 2*Pz_0*z_km1 - 2*y_k*y_km1 - 2*z_k*z_km1 + 2*std::pow(y_km1, 2) + 2*std::pow(z_km1, 2)))/(4*std::sqrt(std::pow(y_k, 2)*std::pow(z_km1, 2) + std::pow(y_km1, 2)*std::pow(z_k, 2) + std::pow(Px_0, 2)*(std::pow(y_k, 2) - 2*y_k*y_km1 + std::pow(y_km1, 2) + std::pow(z_k, 2) - 2*z_k*z_km1 + std::pow(z_km1, 2)) + std::pow(Pz_0, 2)*std::pow(y_k, 2) + std::pow(Pz_0, 2)*std::pow(y_km1, 2) + std::pow(Py_0, 2)*std::pow(z_k, 2) + std::pow(Py_0, 2)*std::pow(z_km1, 2) + std::pow(x_km1, 2)*(std::pow(Py_0, 2) - 2*Py_0*y_k + std::pow(Pz_0, 2) - 2*Pz_0*z_k + std::pow(y_k, 2) + std::pow(z_k, 2)) + std::pow(x_k, 2)*(std::pow(Py_0, 2) - 2*Py_0*y_km1 + std::pow(Pz_0, 2) - 2*Pz_0*z_km1 + std::pow(y_km1, 2) + std::pow(z_km1, 2)) - 2*std::pow(Pz_0, 2)*y_k*y_km1 - 2*Py_0*y_k*std::pow(z_km1, 2) - 2*Py_0*y_km1*std::pow(z_k, 2) - 2*Pz_0*std::pow(y_k, 2)*z_km1 - 2*Pz_0*std::pow(y_km1, 2)*z_k - 2*std::pow(Py_0, 2)*z_k*z_km1 + x_k*x_km1*(2*Py_0*y_k + 2*Py_0*y_km1 + 2*Pz_0*z_k + 2*Pz_0*z_km1 - 2*y_k*y_km1 - 2*z_k*z_km1 - 2*std::pow(Py_0, 2) - 2*std::pow(Pz_0, 2)) + Px_0*x_km1*(2*Py_0*y_k - 2*Py_0*y_km1 + 2*Pz_0*z_k - 2*Pz_0*z_km1 + 2*y_k*y_km1 + 2*z_k*z_km1 - 2*std::pow(y_k, 2) - 2*std::pow(z_k, 2)) - Px_0*x_k*(2*Py_0*y_k - 2*Py_0*y_km1 + 2*Pz_0*z_k - 2*Pz_0*z_km1 - 2*y_k*y_km1 - 2*z_k*z_km1 + 2*std::pow(y_km1, 2) + 2*std::pow(z_km1, 2)) + 2*Pz_0*y_k*y_km1*z_k + 2*Pz_0*y_k*y_km1*z_km1 + 2*Py_0*y_k*z_k*z_km1 + 2*Py_0*y_km1*z_k*z_km1 - 2*y_k*y_km1*z_k*z_km1 - 2*Py_0*Pz_0*y_k*z_k + 2*Py_0*Pz_0*y_k*z_km1 + 2*Py_0*Pz_0*y_km1*z_k - 2*Py_0*Pz_0*y_km1*z_km1)) + (2*x_k*(std::pow(Py_0, 2) - 2*Py_0*y_kp1 + std::pow(Pz_0, 2) - 2*Pz_0*z_kp1 + std::pow(y_kp1, 2) + std::pow(z_kp1, 2)) + x_kp1*(2*Py_0*y_k + 2*Py_0*y_kp1 + 2*Pz_0*z_k + 2*Pz_0*z_kp1 - 2*y_k*y_kp1 - 2*z_k*z_kp1 - 2*std::pow(Py_0, 2) - 2*std::pow(Pz_0, 2)) - Px_0*(2*Py_0*y_k - 2*Py_0*y_kp1 + 2*Pz_0*z_k - 2*Pz_0*z_kp1 - 2*y_k*y_kp1 - 2*z_k*z_kp1 + 2*std::pow(y_kp1, 2) + 2*std::pow(z_kp1, 2)) + (x_kp1*(2*Py_0*y_k - 2*Py_0*y_kp1 + 2*Pz_0*z_k - 2*Pz_0*z_kp1 + 2*y_k*y_kp1 + 2*z_k*z_kp1 - 2*std::pow(y_k, 2) - 2*std::pow(z_k, 2)))/N_p - (x_k*(2*Py_0*y_k - 2*Py_0*y_kp1 + 2*Pz_0*z_k - 2*Pz_0*z_kp1 - 2*y_k*y_kp1 - 2*z_k*z_kp1 + 2*std::pow(y_kp1, 2) + 2*std::pow(z_kp1, 2)))/N_p + (2*Px_0*(std::pow(y_k, 2) - 2*y_k*y_kp1 + std::pow(y_kp1, 2) + std::pow(z_k, 2) - 2*z_k*z_kp1 + std::pow(z_kp1, 2)))/N_p)/(4*std::sqrt(std::pow(y_k, 2)*std::pow(z_kp1, 2) + std::pow(y_kp1, 2)*std::pow(z_k, 2) + std::pow(Px_0, 2)*(std::pow(y_k, 2) - 2*y_k*y_kp1 + std::pow(y_kp1, 2) + std::pow(z_k, 2) - 2*z_k*z_kp1 + std::pow(z_kp1, 2)) + std::pow(Pz_0, 2)*std::pow(y_k, 2) + std::pow(Pz_0, 2)*std::pow(y_kp1, 2) + std::pow(Py_0, 2)*std::pow(z_k, 2) + std::pow(Py_0, 2)*std::pow(z_kp1, 2) + std::pow(x_kp1, 2)*(std::pow(Py_0, 2) - 2*Py_0*y_k + std::pow(Pz_0, 2) - 2*Pz_0*z_k + std::pow(y_k, 2) + std::pow(z_k, 2)) + std::pow(x_k, 2)*(std::pow(Py_0, 2) - 2*Py_0*y_kp1 + std::pow(Pz_0, 2) - 2*Pz_0*z_kp1 + std::pow(y_kp1, 2) + std::pow(z_kp1, 2)) - 2*std::pow(Pz_0, 2)*y_k*y_kp1 - 2*Py_0*y_k*std::pow(z_kp1, 2) - 2*Py_0*y_kp1*std::pow(z_k, 2) - 2*Pz_0*std::pow(y_k, 2)*z_kp1 - 2*Pz_0*std::pow(y_kp1, 2)*z_k - 2*std::pow(Py_0, 2)*z_k*z_kp1 + x_k*x_kp1*(2*Py_0*y_k + 2*Py_0*y_kp1 + 2*Pz_0*z_k + 2*Pz_0*z_kp1 - 2*y_k*y_kp1 - 2*z_k*z_kp1 - 2*std::pow(Py_0, 2) - 2*std::pow(Pz_0, 2)) + Px_0*x_kp1*(2*Py_0*y_k - 2*Py_0*y_kp1 + 2*Pz_0*z_k - 2*Pz_0*z_kp1 + 2*y_k*y_kp1 + 2*z_k*z_kp1 - 2*std::pow(y_k, 2) - 2*std::pow(z_k, 2)) - Px_0*x_k*(2*Py_0*y_k - 2*Py_0*y_kp1 + 2*Pz_0*z_k - 2*Pz_0*z_kp1 - 2*y_k*y_kp1 - 2*z_k*z_kp1 + 2*std::pow(y_kp1, 2) + 2*std::pow(z_kp1, 2)) + 2*Pz_0*y_k*y_kp1*z_k + 2*Pz_0*y_k*y_kp1*z_kp1 + 2*Py_0*y_k*z_k*z_kp1 + 2*Py_0*y_kp1*z_k*z_kp1 - 2*y_k*y_kp1*z_k*z_kp1 - 2*Py_0*Pz_0*y_k*z_k + 2*Py_0*Pz_0*y_k*z_kp1 + 2*Py_0*Pz_0*y_kp1*z_k - 2*Py_0*Pz_0*y_kp1*z_kp1));
    double dAdy = (2*y_k*(std::pow(Px_0, 2) - 2*Px_0*x_km1 + std::pow(Pz_0, 2) - 2*Pz_0*z_km1 + std::pow(x_km1, 2) + std::pow(z_km1, 2)) + y_km1*(2*Px_0*x_k + 2*Px_0*x_km1 + 2*Pz_0*z_k + 2*Pz_0*z_km1 - 2*x_k*x_km1 - 2*z_k*z_km1 - 2*std::pow(Px_0, 2) - 2*std::pow(Pz_0, 2)) - Py_0*(2*Px_0*x_k - 2*Px_0*x_km1 + 2*Pz_0*z_k - 2*Pz_0*z_km1 - 2*x_k*x_km1 - 2*z_k*z_km1 + 2*std::pow(x_km1, 2) + 2*std::pow(z_km1, 2)))/(4*std::sqrt(std::pow(x_k, 2)*std::pow(z_km1, 2) + std::pow(x_km1, 2)*std::pow(z_k, 2) + std::pow(Py_0, 2)*(std::pow(x_k, 2) - 2*x_k*x_km1 + std::pow(x_km1, 2) + std::pow(z_k, 2) - 2*z_k*z_km1 + std::pow(z_km1, 2)) + std::pow(Pz_0, 2)*std::pow(x_k, 2) + std::pow(Pz_0, 2)*std::pow(x_km1, 2) + std::pow(Px_0, 2)*std::pow(z_k, 2) + std::pow(Px_0, 2)*std::pow(z_km1, 2) + std::pow(y_km1, 2)*(std::pow(Px_0, 2) - 2*Px_0*x_k + std::pow(Pz_0, 2) - 2*Pz_0*z_k + std::pow(x_k, 2) + std::pow(z_k, 2)) + std::pow(y_k, 2)*(std::pow(Px_0, 2) - 2*Px_0*x_km1 + std::pow(Pz_0, 2) - 2*Pz_0*z_km1 + std::pow(x_km1, 2) + std::pow(z_km1, 2)) - 2*std::pow(Pz_0, 2)*x_k*x_km1 - 2*Px_0*x_k*std::pow(z_km1, 2) - 2*Px_0*x_km1*std::pow(z_k, 2) - 2*Pz_0*std::pow(x_k, 2)*z_km1 - 2*Pz_0*std::pow(x_km1, 2)*z_k - 2*std::pow(Px_0, 2)*z_k*z_km1 + y_k*y_km1*(2*Px_0*x_k + 2*Px_0*x_km1 + 2*Pz_0*z_k + 2*Pz_0*z_km1 - 2*x_k*x_km1 - 2*z_k*z_km1 - 2*std::pow(Px_0, 2) - 2*std::pow(Pz_0, 2)) + Py_0*y_km1*(2*Px_0*x_k - 2*Px_0*x_km1 + 2*Pz_0*z_k - 2*Pz_0*z_km1 + 2*x_k*x_km1 + 2*z_k*z_km1 - 2*std::pow(x_k, 2) - 2*std::pow(z_k, 2)) - Py_0*y_k*(2*Px_0*x_k - 2*Px_0*x_km1 + 2*Pz_0*z_k - 2*Pz_0*z_km1 - 2*x_k*x_km1 - 2*z_k*z_km1 + 2*std::pow(x_km1, 2) + 2*std::pow(z_km1, 2)) + 2*Pz_0*x_k*x_km1*z_k + 2*Pz_0*x_k*x_km1*z_km1 + 2*Px_0*x_k*z_k*z_km1 + 2*Px_0*x_km1*z_k*z_km1 - 2*x_k*x_km1*z_k*z_km1 - 2*Px_0*Pz_0*x_k*z_k + 2*Px_0*Pz_0*x_k*z_km1 + 2*Px_0*Pz_0*x_km1*z_k - 2*Px_0*Pz_0*x_km1*z_km1)) + (2*y_k*(std::pow(Px_0, 2) - 2*Px_0*x_kp1 + std::pow(Pz_0, 2) - 2*Pz_0*z_kp1 + std::pow(x_kp1, 2) + std::pow(z_kp1, 2)) + y_kp1*(2*Px_0*x_k + 2*Px_0*x_kp1 + 2*Pz_0*z_k + 2*Pz_0*z_kp1 - 2*x_k*x_kp1 - 2*z_k*z_kp1 - 2*std::pow(Px_0, 2) - 2*std::pow(Pz_0, 2)) - Py_0*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Pz_0*z_k - 2*Pz_0*z_kp1 - 2*x_k*x_kp1 - 2*z_k*z_kp1 + 2*std::pow(x_kp1, 2) + 2*std::pow(z_kp1, 2)) + (y_kp1*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Pz_0*z_k - 2*Pz_0*z_kp1 + 2*x_k*x_kp1 + 2*z_k*z_kp1 - 2*std::pow(x_k, 2) - 2*std::pow(z_k, 2)))/N_p - (y_k*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Pz_0*z_k - 2*Pz_0*z_kp1 - 2*x_k*x_kp1 - 2*z_k*z_kp1 + 2*std::pow(x_kp1, 2) + 2*std::pow(z_kp1, 2)))/N_p + (2*Py_0*(std::pow(x_k, 2) - 2*x_k*x_kp1 + std::pow(x_kp1, 2) + std::pow(z_k, 2) - 2*z_k*z_kp1 + std::pow(z_kp1, 2)))/N_p)/(4*std::sqrt(std::pow(x_k, 2)*std::pow(z_kp1, 2) + std::pow(x_kp1, 2)*std::pow(z_k, 2) + std::pow(Py_0, 2)*(std::pow(x_k, 2) - 2*x_k*x_kp1 + std::pow(x_kp1, 2) + std::pow(z_k, 2) - 2*z_k*z_kp1 + std::pow(z_kp1, 2)) + std::pow(Pz_0, 2)*std::pow(x_k, 2) + std::pow(Pz_0, 2)*std::pow(x_kp1, 2) + std::pow(Px_0, 2)*std::pow(z_k, 2) + std::pow(Px_0, 2)*std::pow(z_kp1, 2) + std::pow(y_kp1, 2)*(std::pow(Px_0, 2) - 2*Px_0*x_k + std::pow(Pz_0, 2) - 2*Pz_0*z_k + std::pow(x_k, 2) + std::pow(z_k, 2)) + std::pow(y_k, 2)*(std::pow(Px_0, 2) - 2*Px_0*x_kp1 + std::pow(Pz_0, 2) - 2*Pz_0*z_kp1 + std::pow(x_kp1, 2) + std::pow(z_kp1, 2)) - 2*std::pow(Pz_0, 2)*x_k*x_kp1 - 2*Px_0*x_k*std::pow(z_kp1, 2) - 2*Px_0*x_kp1*std::pow(z_k, 2) - 2*Pz_0*std::pow(x_k, 2)*z_kp1 - 2*Pz_0*std::pow(x_kp1, 2)*z_k - 2*std::pow(Px_0, 2)*z_k*z_kp1 + y_k*y_kp1*(2*Px_0*x_k + 2*Px_0*x_kp1 + 2*Pz_0*z_k + 2*Pz_0*z_kp1 - 2*x_k*x_kp1 - 2*z_k*z_kp1 - 2*std::pow(Px_0, 2) - 2*std::pow(Pz_0, 2)) + Py_0*y_kp1*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Pz_0*z_k - 2*Pz_0*z_kp1 + 2*x_k*x_kp1 + 2*z_k*z_kp1 - 2*std::pow(x_k, 2) - 2*std::pow(z_k, 2)) - Py_0*y_k*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Pz_0*z_k - 2*Pz_0*z_kp1 - 2*x_k*x_kp1 - 2*z_k*z_kp1 + 2*std::pow(x_kp1, 2) + 2*std::pow(z_kp1, 2)) + 2*Pz_0*x_k*x_kp1*z_k + 2*Pz_0*x_k*x_kp1*z_kp1 + 2*Px_0*x_k*z_k*z_kp1 + 2*Px_0*x_kp1*z_k*z_kp1 - 2*x_k*x_kp1*z_k*z_kp1 - 2*Px_0*Pz_0*x_k*z_k + 2*Px_0*Pz_0*x_k*z_kp1 + 2*Px_0*Pz_0*x_kp1*z_k - 2*Px_0*Pz_0*x_kp1*z_kp1));
    double dAdz = (2*z_k*(std::pow(Px_0, 2) - 2*Px_0*x_kp1 + std::pow(Py_0, 2) - 2*Py_0*y_kp1 + std::pow(x_kp1, 2) + std::pow(y_kp1, 2)) + z_kp1*(2*Px_0*x_k + 2*Px_0*x_kp1 + 2*Py_0*y_k + 2*Py_0*y_kp1 - 2*x_k*x_kp1 - 2*y_k*y_kp1 - 2*std::pow(Px_0, 2) - 2*std::pow(Py_0, 2)) - Pz_0*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Py_0*y_k - 2*Py_0*y_kp1 - 2*x_k*x_kp1 - 2*y_k*y_kp1 + 2*std::pow(x_kp1, 2) + 2*std::pow(y_kp1, 2)) + (z_kp1*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Py_0*y_k - 2*Py_0*y_kp1 + 2*x_k*x_kp1 + 2*y_k*y_kp1 - 2*std::pow(x_k, 2) - 2*std::pow(y_k, 2)))/N_p - (z_k*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Py_0*y_k - 2*Py_0*y_kp1 - 2*x_k*x_kp1 - 2*y_k*y_kp1 + 2*std::pow(x_kp1, 2) + 2*std::pow(y_kp1, 2)))/N_p + (2*Pz_0*(std::pow(x_k, 2) - 2*x_k*x_kp1 + std::pow(x_kp1, 2) + std::pow(y_k, 2) - 2*y_k*y_kp1 + std::pow(y_kp1, 2)))/N_p)/(4*std::sqrt(std::pow(x_k, 2)*std::pow(y_kp1, 2) + std::pow(x_kp1, 2)*std::pow(y_k, 2) + std::pow(Pz_0, 2)*(std::pow(x_k, 2) - 2*x_k*x_kp1 + std::pow(x_kp1, 2) + std::pow(y_k, 2) - 2*y_k*y_kp1 + std::pow(y_kp1, 2)) + std::pow(Py_0, 2)*std::pow(x_k, 2) + std::pow(Py_0, 2)*std::pow(x_kp1, 2) + std::pow(Px_0, 2)*std::pow(y_k, 2) + std::pow(Px_0, 2)*std::pow(y_kp1, 2) + std::pow(z_kp1, 2)*(std::pow(Px_0, 2) - 2*Px_0*x_k + std::pow(Py_0, 2) - 2*Py_0*y_k + std::pow(x_k, 2) + std::pow(y_k, 2)) + std::pow(z_k, 2)*(std::pow(Px_0, 2) - 2*Px_0*x_kp1 + std::pow(Py_0, 2) - 2*Py_0*y_kp1 + std::pow(x_kp1, 2) + std::pow(y_kp1, 2)) - 2*std::pow(Py_0, 2)*x_k*x_kp1 - 2*Px_0*x_k*std::pow(y_kp1, 2) - 2*Px_0*x_kp1*std::pow(y_k, 2) - 2*Py_0*std::pow(x_k, 2)*y_kp1 - 2*Py_0*std::pow(x_kp1, 2)*y_k - 2*std::pow(Px_0, 2)*y_k*y_kp1 + z_k*z_kp1*(2*Px_0*x_k + 2*Px_0*x_kp1 + 2*Py_0*y_k + 2*Py_0*y_kp1 - 2*x_k*x_kp1 - 2*y_k*y_kp1 - 2*std::pow(Px_0, 2) - 2*std::pow(Py_0, 2)) + Pz_0*z_kp1*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Py_0*y_k - 2*Py_0*y_kp1 + 2*x_k*x_kp1 + 2*y_k*y_kp1 - 2*std::pow(x_k, 2) - 2*std::pow(y_k, 2)) - Pz_0*z_k*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Py_0*y_k - 2*Py_0*y_kp1 - 2*x_k*x_kp1 - 2*y_k*y_kp1 + 2*std::pow(x_kp1, 2) + 2*std::pow(y_kp1, 2)) + 2*Py_0*x_k*x_kp1*y_k + 2*Py_0*x_k*x_kp1*y_kp1 + 2*Px_0*x_k*y_k*y_kp1 + 2*Px_0*x_kp1*y_k*y_kp1 - 2*x_k*x_kp1*y_k*y_kp1 - 2*Px_0*Py_0*x_k*y_k + 2*Px_0*Py_0*x_k*y_kp1 + 2*Px_0*Py_0*x_kp1*y_k - 2*Px_0*Py_0*x_kp1*y_kp1)) + (2*z_k*(std::pow(Px_0, 2) - 2*Px_0*x_km1 + std::pow(Py_0, 2) - 2*Py_0*y_km1 + std::pow(x_km1, 2) + std::pow(y_km1, 2)) + z_km1*(2*Px_0*x_k + 2*Px_0*x_km1 + 2*Py_0*y_k + 2*Py_0*y_km1 - 2*x_k*x_km1 - 2*y_k*y_km1 - 2*std::pow(Px_0, 2) - 2*std::pow(Py_0, 2)) - Pz_0*(2*Px_0*x_k - 2*Px_0*x_km1 + 2*Py_0*y_k - 2*Py_0*y_km1 - 2*x_k*x_km1 - 2*y_k*y_km1 + 2*std::pow(x_km1, 2) + 2*std::pow(y_km1, 2)))/(4*std::sqrt(std::pow(x_k, 2)*std::pow(y_km1, 2) + std::pow(x_km1, 2)*std::pow(y_k, 2) + std::pow(Pz_0, 2)*(std::pow(x_k, 2) - 2*x_k*x_km1 + std::pow(x_km1, 2) + std::pow(y_k, 2) - 2*y_k*y_km1 + std::pow(y_km1, 2)) + std::pow(Py_0, 2)*std::pow(x_k, 2) + std::pow(Py_0, 2)*std::pow(x_km1, 2) + std::pow(Px_0, 2)*std::pow(y_k, 2) + std::pow(Px_0, 2)*std::pow(y_km1, 2) + std::pow(z_km1, 2)*(std::pow(Px_0, 2) - 2*Px_0*x_k + std::pow(Py_0, 2) - 2*Py_0*y_k + std::pow(x_k, 2) + std::pow(y_k, 2)) + std::pow(z_k, 2)*(std::pow(Px_0, 2) - 2*Px_0*x_km1 + std::pow(Py_0, 2) - 2*Py_0*y_km1 + std::pow(x_km1, 2) + std::pow(y_km1, 2)) - 2*std::pow(Py_0, 2)*x_k*x_km1 - 2*Px_0*x_k*std::pow(y_km1, 2) - 2*Px_0*x_km1*std::pow(y_k, 2) - 2*Py_0*std::pow(x_k, 2)*y_km1 - 2*Py_0*std::pow(x_km1, 2)*y_k - 2*std::pow(Px_0, 2)*y_k*y_km1 + z_k*z_km1*(2*Px_0*x_k + 2*Px_0*x_km1 + 2*Py_0*y_k + 2*Py_0*y_km1 - 2*x_k*x_km1 - 2*y_k*y_km1 - 2*std::pow(Px_0, 2) - 2*std::pow(Py_0, 2)) + Pz_0*z_km1*(2*Px_0*x_k - 2*Px_0*x_km1 + 2*Py_0*y_k - 2*Py_0*y_km1 + 2*x_k*x_km1 + 2*y_k*y_km1 - 2*std::pow(x_k, 2) - 2*std::pow(y_k, 2)) - Pz_0*z_k*(2*Px_0*x_k - 2*Px_0*x_km1 + 2*Py_0*y_k - 2*Py_0*y_km1 - 2*x_k*x_km1 - 2*y_k*y_km1 + 2*std::pow(x_km1, 2) + 2*std::pow(y_km1, 2)) + 2*Py_0*x_k*x_km1*y_k + 2*Py_0*x_k*x_km1*y_km1 + 2*Px_0*x_k*y_k*y_km1 + 2*Px_0*x_km1*y_k*y_km1 - 2*x_k*x_km1*y_k*y_km1 - 2*Px_0*Py_0*x_k*y_k + 2*Px_0*Py_0*x_k*y_km1 + 2*Px_0*Py_0*x_km1*y_k - 2*Px_0*Py_0*x_km1*y_km1));
    std::array<double, 3> dAdr = {dAdx, dAdy, dAdz};
    return dAdr; 
}
std::array<double, 3> Simulation::dVdr(Vertex* current, Vertex* prev, Vertex* next, std::array<double,3> polyCenter, std::array<double,3> cellCenter, int N_p, int N_c){
    double x_k = current->getPos()[0];
    double y_k = current->getPos()[1];
    double z_k = current->getPos()[2];

    double x_kp1 = next->getPos()[0];
    double y_kp1 = next->getPos()[1];
    double z_kp1 = next->getPos()[2];

    double x_km1 = prev->getPos()[0];
    double y_km1 = prev->getPos()[1];
    double z_km1 = prev->getPos()[2];

    double Px_0 = polyCenter[0];
    double Py_0 = polyCenter[1];
    double Pz_0 = polyCenter[2];

    double Cx_0 = cellCenter[0];
    double Cy_0 = cellCenter[1];
    double Cz_0 = cellCenter[2]; 
    
    double dVdx = ((Cz_0 - Pz_0)*(Py_0 - y_kp1))/6 - ((Cz_0 - Pz_0)*(Py_0 - y_km1))/6 - ((Cz_0 - Pz_0)*(y_k - y_kp1) - (Cy_0 - Py_0)*(z_k - z_kp1) - (Py_0 - y_k)*(Pz_0 - z_kp1) + (Py_0 - y_kp1)*(Pz_0 - z_k))/(6*N_p) + ((Cy_0 - Py_0)*(Pz_0 - z_km1))/6 - ((Cy_0 - Py_0)*(Pz_0 - z_kp1))/6 - ((Py_0 - y_k)*(Pz_0 - z_kp1) - (Py_0 - y_kp1)*(Pz_0 - z_k))/(6*N_c);
    double dVdy = ((Cz_0 - Pz_0)*(x_k - x_kp1) - (Cx_0 - Px_0)*(z_k - z_kp1) - (Px_0 - x_k)*(Pz_0 - z_kp1) + (Px_0 - x_kp1)*(Pz_0 - z_k))/(6*N_p) + ((Cz_0 - Pz_0)*(Px_0 - x_km1))/6 - ((Cz_0 - Pz_0)*(Px_0 - x_kp1))/6 - ((Cx_0 - Px_0)*(Pz_0 - z_km1))/6 + ((Cx_0 - Px_0)*(Pz_0 - z_kp1))/6 + ((Px_0 - x_k)*(Pz_0 - z_kp1) - (Px_0 - x_kp1)*(Pz_0 - z_k))/(6*N_c);
    double dVdz = ((Cy_0 - Py_0)*(Px_0 - x_kp1))/6 - ((Cy_0 - Py_0)*(Px_0 - x_km1))/6 - ((Cy_0 - Py_0)*(x_k - x_kp1) - (Cx_0 - Px_0)*(y_k - y_kp1) - (Px_0 - x_k)*(Py_0 - y_kp1) + (Px_0 - x_kp1)*(Py_0 - y_k))/(6*N_p) + ((Cx_0 - Px_0)*(Py_0 - y_km1))/6 - ((Cx_0 - Px_0)*(Py_0 - y_kp1))/6 - ((Px_0 - x_k)*(Py_0 - y_kp1) - (Px_0 - x_kp1)*(Py_0 - y_k))/(6*N_c);

    std::array<double, 3> dVdr = {dVdx, dVdy, dVdz};
    return dVdr;
}

std::string Simulation::convertDouble(double val){
    std::string str = std::to_string(val);  // Convert double to string
    
    // Remove trailing zeros
    str.erase(str.find_last_not_of('0') + 1, std::string::npos);
    
    // Remove decimal point if it's the last character
    if (str.back() == '.') {
        str.pop_back();
    }
    std::replace(str.begin(), str.end(), '.', 'p');
    return str;
}