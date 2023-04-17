#include "system.h"
using namespace Eigen;

float System::updateWaterGrid() {
    float timeStep = calcTimeStep();
    updateGridFromMarkers();
    updateVelocityField();
    return timeStep;
}

float System::calcTimeStep() {
    return 0;
}


void System::updateGridFromMarkers() {
    // TODO
}



void System::updateVelocityField() {
    // TODO
}
