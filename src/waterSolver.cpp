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

}



void System::updateVelocityField() {

}
