#include "system.h"
using namespace Eigen;
using namespace std;

void System::solve(){
    float timeStep = updateWaterGrid();
    updateParticles(timeStep);
}

void System::updateParticles(float timeStep){
    checkNanAndInf();
    #pragma omp parallel for
    for (Particle &inkPtcl: m_ink) {
        // equation 9
        inkPtcl.velocity = DENSITY * getVelocity(inkPtcl.position);
        // equation 7
        Vector3f midParticlePos = inkPtcl.position + (inkPtcl.velocity * timeStep * .5);
        // get midpoint particle pos from velocity field
        Vector3f midParticleVel = DENSITY * getVelocity(midParticlePos);
        // equation 8
        inkPtcl.position = midParticlePos + (timeStep * midParticleVel);

        // Fixing x
        if (inkPtcl.position.x() < 0 || inkPtcl.position.x() > WATERGRID_X*CELL_DIM) {
            inkPtcl.velocity[0] *= -1.f;
            inkPtcl.position[0] = min(max(inkPtcl.position[0], 0.01f), WATERGRID_X*CELL_DIM-0.01f);
        }

        // Fixing y
        if (inkPtcl.position.y() < 0 || inkPtcl.position.y() > WATERGRID_Y*CELL_DIM) {
            inkPtcl.velocity[1] *= -1.f;
            inkPtcl.position[1] = min(max(inkPtcl.position[1], 0.01f), WATERGRID_Y*CELL_DIM-0.01f);
        }

        // Fixing y
        if (inkPtcl.position.z() < 0 || inkPtcl.position.z() > WATERGRID_Z*CELL_DIM) {
            inkPtcl.velocity[2] *= -1.f;
            inkPtcl.position[2] = min(max(inkPtcl.position[2], 0.01f), WATERGRID_Z*CELL_DIM-0.01f);
        }


        // opacity
        // lifetime
    }
}
