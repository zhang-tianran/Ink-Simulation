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

Vector3i System::getCellIndexFromPoint(Vector3f &pos) {
    return Vector3i{floor(pos.x()), floor(pos.y()), floor(pos.z())};
}

/// Applies the convection term in the Navier-Stokes equation to each cell's velocity
void System::applyConvection(float timeStep) {
    for (auto &kv : m_waterGrid) {
        assert(kv.second.oldVelocity == kv.second.currVelocity);
    }

    /// Update each cell's curr_velocity in the waterGrid...
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                Cell currCell = m_waterGrid.at(Vector3i{i, j, k});
                Vector3f virtualParticlePos = traceParticle(i + CELL_DIM/2.f, j + CELL_DIM/2.f, k + CELL_DIM/2.f, -timeStep);
                currCell.currVelocity = m_waterGrid.at(getCellIndexFromPoint(virtualParticlePos)).oldVelocity;
            }
        }
    }

    /// Update each cell's old_velocity to be the curr_velocity
    for (auto &kv : m_waterGrid) {
        kv.second.oldVelocity = kv.second.currVelocity;
    }
}

/// Applies the external force term in the Navier-Stokes equation to each cell's velocity
void System::applyExternalForces(float timeStep) {
    for (auto &kv : m_waterGrid) {
        /// Applies gravity
        kv.second.currVelocity += gravity * timeStep;
        kv.second.oldVelocity = kv.second.currVelocity;

        // TODO: Add vorticity confinement force
    }
}

/**
 * Performs the Laplacian Operator on the velocity vector field.
 *
 * i, j, k : cell index
 * idx     : component of the velocity vector (either 0, 1, or 2 for x, y, or z)
 */
float System::laplacianOperatorOnVelocity(int i, int j, int k, int idx) {
    float laplacianVelocity = 0;

    /// i direction
    laplacianVelocity += (i+1 < WATERGRID_X) ? m_waterGrid.at(Vector3i{i+1, j, k}).oldVelocity[idx] : 0;
    laplacianVelocity += (i-1 >= 0         ) ? m_waterGrid.at(Vector3i{i-1, j, k}).oldVelocity[idx] : 0;

    /// j direction
    laplacianVelocity += (j+1 < WATERGRID_Y) ? m_waterGrid.at(Vector3i{i, j+1, k}).oldVelocity[idx] : 0;
    laplacianVelocity += (j-1 >= 0         ) ? m_waterGrid.at(Vector3i{i, j-1, k}).oldVelocity[idx] : 0;

    /// k direction
    laplacianVelocity += (k+1 < WATERGRID_Z) ? m_waterGrid.at(Vector3i{i, j, k+1}).oldVelocity[idx] : 0;
    laplacianVelocity += (k-1 >= 0         ) ? m_waterGrid.at(Vector3i{i, j, k-1}).oldVelocity[idx] : 0;

    /// -6*currCellOldVelocity term
    laplacianVelocity -= 6 * m_waterGrid.at(Vector3i{i, j, k}).oldVelocity[idx];

    return laplacianVelocity;
}

/// Applies the viscosity term in the Navier-Stokes equation to each cell's velocity
void System::applyViscosity(float timeStep) {
    /// For each cell, apply the viscosity to each cell's currVelocity
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                float u_x = laplacianOperatorOnVelocity(i, j, k, 0);
                float u_y = laplacianOperatorOnVelocity(i, j, k, 1);
                float u_z = laplacianOperatorOnVelocity(i, j, k, 2);
                m_waterGrid.at(Vector3f{i, j, k}).currVelocity += Vector3f{u_x, u_y, u_z};
            }
        }
    }

    /// Make the oldVelocity = currVelocity
    for (auto &kv : m_waterGrid) {
        kv.second.oldVelocity = kv.second.currVelocity;
    }
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

/// AP = B (equation 13)
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

/// Applies the pressure term in the Navier-Stokes equation to each cell's velocity
void System::applyPressure(float timeStep) {
    MatrixXf pressure = calculatePressure(timeStep);
    for (int l = 0; l < WATERGRID_X; l++) {
        for (int w = 0; w < WATERGRID_Z; w++) {
            for (int h = 0; h < WATERGRID_Y; h++) {
                int row_idx = grid2mat(l, w, h);
                Vector3f gradient = getGradient(l, w, h, pressure);
                m_waterGrid[Vector3i(l, w, h)].currVelocity -= (timeStep / DENSITY * CELL_DIM) * gradient;
            }
        }
    }
}
