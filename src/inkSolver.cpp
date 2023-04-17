#include "system.h"
using namespace Eigen;

void System::solve(){
    float timeStep = updateWaterGrid();
    updateParticles(timeStep);
}

void System::updateParticles(float timeStep){
    for (Particle inkPtcl: m_ink) {
        // midpoint update
        Vector3f midCellVel; // = TODO: cell velocity of grid correspnding to current position
        Vector3f midPos = inkPtcl.position + inkPtcl.velocity * timeStep / 2.f;
        Vector3f midVel = midCellVel * timeStep / 2.f;
        // final update
        inkPtcl.position = midPos + timeStep * midVel;
        inkPtcl.velocity = midCellVel * timeStep;

        // opacity
        // lifetime
    }
}
