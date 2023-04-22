#include "system.h"
using namespace Eigen;
using namespace std;

void System::solve(){
    float timeStep = updateWaterGrid();
    updateParticles(timeStep);
}

void System::updateParticles(float timeStep){
    checkNanAndInf();
    for (Particle &inkPtcl: m_ink) {
//        // midpoint velocity
//        Vector3f midVel = getVelocity(inkPtcl.position);
//        Vector3f midPos = inkPtcl.position + (midVel * timeStep * .5);
//        Vector3f finalVel = getVelocity(midPos);
//        // final update
//        inkPtcl.position += timeStep * finalVel;
//        inkPtcl.velocity = finalVel; // TODO: this field is never used...

        // equation 7
        Vector3f midParticlePos = inkPtcl.position + (inkPtcl.velocity * timeStep * .5);
        // get midpoint particle pos from velocity field
        Vector3f midParticleVel = DENSITY*getVelocity(midParticlePos);
        // equation 9
        inkPtcl.velocity = DENSITY*getVelocity(inkPtcl.position);
        // equation 8
        inkPtcl.position = midParticlePos + (timeStep*midParticleVel);

        // opacity
        // lifetime
    }

    int i = 1;
}
