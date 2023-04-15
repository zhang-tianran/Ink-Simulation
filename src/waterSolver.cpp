#include "system.h"
using namespace Eigen;

float System::updateWaterGrid(std::unordered_map<Eigen::Vector3f, Cell, hash_func> &waterGrid) {
    float timeStep = calcTimeStep(waterGrid);
    updateGridFromMarkers(waterGrid);
    updateVelocityField(waterGrid);
    return timeStep;
}

float System::calcTimeStep(const std::unordered_map<Eigen::Vector3f, Cell, hash_func> &waterGrid) {
    return 0;
}


void System::updateGridFromMarkers(std::unordered_map<Eigen::Vector3f, Cell, hash_func> &waterGrid) {

}



void System::updateVelocityField(std::unordered_map<Eigen::Vector3f, Cell, hash_func> &waterGrid) {

}
