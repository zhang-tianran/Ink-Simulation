#include "system.h"
using namespace Eigen;

void System::solve(){
    float timeStep = updateWaterGrid();
    updateParticles(timeStep);
}

void System::updateParticles(float timeStep){
    for (Particle &inkPtcl: m_ink) {
        // midpoint velocity
        Vector3f midVel = getVelocity(inkPtcl.position);
        Vector3f midPos = inkPtcl.position + midVel * timeStep / 2.f;
        Vector3f finalVel = getVelocity(midPos);
        // final update
        inkPtcl.position += timeStep * finalVel;
        inkPtcl.velocity = finalVel;

        // opacity
        // lifetime
    }
}
