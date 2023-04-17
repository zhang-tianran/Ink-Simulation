#include "system.h"
using namespace Eigen;

float System::updateWaterGrid() {
    float timeStep = calcTimeStep();
    updateVelocityField(timeStep);
    return timeStep;
}

float System::calcTimeStep() {
    return 0;
}

/// See "Fluid Flow 4 the Rest of Us" Paper for more details
void System::updateVelocityField(float timeStep) {
    applyConvection(timeStep);
    applyExternalForces(timeStep);
    applyViscosity(timeStep);
    calculatePressure();
    applyPressure(timeStep);
}

void System::applyConvection(float timeStep) {
    // TODO
}

void System::applyExternalForces(float timeStep) {
    // TODO
}

void System::applyViscosity(float timeStep) {
    // TODO
}

void System::calculatePressure() {
    // TODO
}

void System::applyPressure(float timeStep) {
    // TODO
}
