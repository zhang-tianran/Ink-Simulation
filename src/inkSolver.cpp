#include "system.h"
using namespace Eigen;
using namespace std;

double System::solve(double timeToNextRender){
    double timeStep = updateWaterGrid(timeToNextRender);
    updateParticles(timeStep);
    return timeStep;
}

void System::emitParticleHemisphere(float radius) {
    // emit outwards from a hemisphere
    float fraction = M_PI / 12.f;
    float twoPi = M_PI * 2.f;
    Vector3f surfaceNormal = Vector3f(WATERGRID_X / 2.f, -WATERGRID_Y -.2f, WATERGRID_Z / 2.f); // sampling point on middle of ceiling, pointing downwards
    Vector3f up = Vector3f(0, 1, 0);
    Vector3f axis = (surfaceNormal + up).normalized();
    AngleAxisf rot(M_PI, axis);
    for (float i = 0; i < twoPi; i += fraction) {
        for (float j = 0; j < M_PI; j += fraction) {
            float x = radius*sin(j) * sin(i);
            float y = radius*cos(j); //
            float z = radius*sin(j) * cos(i);

            // starting position of particle
            Vector3f pos = Vector3f(x, y, z);
            // starting velocity of particle
            Vector3f vel = rot * pos;
            Particle p {
                .position = pos,
                .velocity = vel,
                .opacity = 1.f,
                .lifeTime = 5.f,
            };

            this->m_ink.push_back(p);
        }
    }

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


        if (USE_LIFETIME) {
            inkPtcl.lifeTime -= timeStep;
            // TODO: idk if this is good, but it assumes starting lifetime is 5f
            inkPtcl.opacity = inkPtcl.lifeTime / 5.f;
        }
        // opacity
        // lifetime
    }
}
