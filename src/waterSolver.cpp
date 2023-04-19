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
    applyPressure(timeStep);
}

Vector3f System::traceParticle(float x, float y, float z, float t) {
    Vector3f vel = getVelocity(Vector3f(x, y, z));
    vel = getVelocity(Vector3f(x + .5*t*vel.x(), y + .5*t*vel.y(), z + .5*t*vel.z()));
    return Vector3f(x, y, z) + t*vel;
}

Vector3i System::getCellIndexFromPoint(float x, float y, float z) {
    return Vector3i{floor(x), floor(y), floor(z)};
}

//// Applies the convection term in the Navier-Stokes equation
void System::applyConvection(float timeStep) {
    for (auto &kv : m_waterGrid) {
        assert(kv.second.old_velocity == kv.second.curr_velocity);
    }

    ///// Update each cell's curr_velocity in the waterGrid...
    for (int x = 0; x < WATERGRID_X; x++) {
        for (int y = 0; y < WATERGRID_Y; y++) {
            for (int z = 0; z < WATERGRID_Z; z++) {
                Cell currCell = m_waterGrid.at(Vector3i{x, y, z});
                Vector3f virtualParticlePos = traceParticle(x + CELL_DIM/2.f, y + CELL_DIM/2.f, z + CELL_DIM/2.f);
                currCell.curr_velocity = m_waterGrid.at(getCellIndexFromPoint(virtualParticlePos)).old_velocity;
            }
        }
    }

    //// Update each cell's old_velocity to be the curr_velocity
    for (auto &kv : m_waterGrid) {
        kv.second.old_velocity = kv.second.curr_velocity;
    }
}

//// Applies the external force term in the Navier-Stokes equation
void System::applyExternalForces(float timeStep) {
    for (auto& kv : m_waterGrid) {
        //// Applies gravity
        kv.second.curr_velocity = kv.second.curr_velocity + gravity * timeStep; // TODO: double check if we can just add gravity like this
        kv.second.old_velocity = kv.second.curr_velocity;

        // TODO: Add vorticity confinement force
    }
}

void System::applyViscosity(float timeStep) {
    // TODO
    // for each cell, see equation 7 and 6
    // for u_x, we get the x component of neighboring cells and do equation 6
}

void System::initPressureA() {
    int n = WATERGRID_X * WATERGRID_Y * WATERGRID_Z;
    SpMat A(n, n);
    A.set_zeros();
    for (int l = 0; l < WATERGRID_X; l++) {
        for (int w = 0; w < WATERGRID_Z; w++) {
            for (int h = 0; h < WATERGRID_Y; h++) {
                int row_idx = grid2mat(l, w, h);
                std::vector<Vector3i> neighbors = getGridNeighbors(l, w, h);
                for (Vector3i neighbor : neighbors) {
                    A.insert(row_idx, grid2mat(neighbor[0], neighbor[1], neighbor[2]));
                }

                A.insert(row_idx, row_idx) = -6;
            }
        }
    }
    llt.compute(A);
}

// AP = B (equation 13)
MatrixXf System::calculatePressure(float timeStep) {
    int n = WATERGRID_X * WATERGRID_Y * WATERGRID_Z;
    MatrixXf bMatrix(n, 1);
    for (int l = 0; l < WATERGRID_X; l++) {
        for (int w = 0; w < WATERGRID_Z; w++) {
            for (int h = 0; h < WATERGRID_Y; h++) {
                int row_idx = grid2mat(l, w, h);
                float divergence =  getDivergence(l, w, h);
                std::vector<Vector3i> neighbors = getGridNeighbors(l, w, h);
                int k = 6 - neighbors.size();
                bMatrix.row(row_idx) = ((DENSITY * CELL_DIM) / timeStep) * divergence - k * ATMOSPHERIC_PRESSURE; 
            }
        }
    }
    MatrixXf pressure = llt.solve(bMatrix);
    return pressure;
}

void System::applyPressure(float timeStep) {
    MatrixXf pressure = calculatePressure(timeStep);
    for (int l = 0; l < WATERGRID_X; l++) {
        for (int w = 0; w < WATERGRID_Z; w++) {
            for (int h = 0; h < WATERGRID_Y; h++) {
                int row_idx = grid2mat(l, w, h);
                Vector3f gradient = getGradient(l, w, h, pressure);
                m_waterGrid[Vector3i(l, w, h)].velocity -= (timeStep / DENSITY * CELL_DIM) * gradient;
            }
        }
    }
}
