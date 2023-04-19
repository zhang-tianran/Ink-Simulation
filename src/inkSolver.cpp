#include "system.h"
using namespace Eigen;

void System::solve(){
    float timeStep = updateWaterGrid();
    updateParticles(timeStep);
}

void System::updateParticles(float timeStep){
    for (Particle inkPtcl: m_ink) {
        // midpoint velocity
        Vector3f midVel = getVelocity(inkPtcl.position);
        Vector3f midPos = inkPtcl.position + inkPtcl.velocity * timeStep / 2.f;
        Vector3f finalVel = getVelocity(midPos);
        // final update
        inkPtcl.position += timeStep * finalVel;
        inkPtcl.velocity = finalVel;

        // opacity
        // lifetime
    }
}

//// Get the interpolated velocity at a point in space.
Vector3f System::getVelocity(Vector3f pos){
    float x = getInterpolatedValue(pos[0] / CELL_DIM, pos[1] / CELL_DIM - 0.5f, pos[2] / CELL_DIM - 0.5f, 0);
    float y = getInterpolatedValue(pos[0] / CELL_DIM - 0.5f, pos[1] / CELL_DIM, pos[2] / CELL_DIM - 0.5f, 1);
    float z = getInterpolatedValue(pos[0] / CELL_DIM - 0.5f, pos[1] / CELL_DIM - 0.5f, pos[2] / CELL_DIM, 2);
    return Vector3f(x, y, z);
}

//// Get an interpolated data value from the grid
float System::getInterpolatedValue(float x, float y, float z, int idx) {
    int i = floor(x);
    int j = floor(y);
    int k = floor(z);
    return (i + 1 - x) * (j + 1 - y) * (k + 1 - z) * m_waterGrid[Vector3i(i, j, k)].velocity[idx] +
           (x - i) * (j + 1 - y) * (k + 1 - z) * m_waterGrid[Vector3i(i + 1, j, k)].velocity[idx] +
           (i + 1 - x) * (y - j) * (k + 1 - z) * m_waterGrid[Vector3i(i, j + 1, k)].velocity[idx] +
           (x - i) * (y - j) * (k + 1 - z) * m_waterGrid[Vector3i(i + 1, j + 1, k)].velocity[idx] +
           (i + 1 - x) * (j + 1 - y) * (z - k) * m_waterGrid[Vector3i(i, j, k + 1)].velocity[idx] +
           (x - i) * (j + 1 - y) * (z - k) * m_waterGrid[Vector3i(i + 1, j, k + 1)].velocity[idx] +
           (i + 1 - x) * (y - j) * (z - k) * m_waterGrid[Vector3i(i, j + 1, k + 1)].velocity[idx] +
           (x - i) * (y - k) * (z - k) * m_waterGrid[Vector3i(i + 1, j + 1, k + 1)].velocity[idx];
}

