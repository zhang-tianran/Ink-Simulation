#include "system.h"

using namespace Eigen;

float System::updateWaterGrid() {
    float timeStep = calcTimeStep();
    updateVelocityField(timeStep);
    return timeStep;
}

float System::calcTimeStep() {
    float maxVelocity = 0;
    for (auto kv : m_waterGrid) {
        if (maxVelocity < kv.second.oldVelocity.norm()) {
            maxVelocity = kv.second.oldVelocity.norm();
        }
    }
    if (maxVelocity == 0) {
        return MIN_TIMESTEP;
    } else {
        float timeStep = K_CFL * (CELL_DIM / maxVelocity);;
        timeStep = std::max(std::min(timeStep, MAX_TIMESTEP), MIN_TIMESTEP);
        return timeStep;
    }
}

/// See "Fluid Flow 4 the Rest of Us" Paper for more details
void System::updateVelocityField(float timeStep) {
    /// Error Checking
    for (auto &kv : m_waterGrid) {
        assert(kv.second.oldVelocity == kv.second.currVelocity);
    }

    /// Navier-Stokes equation
    applyConvection(timeStep);
    applyExternalForces(timeStep);
    applyViscosity(timeStep);
    applyPressure(timeStep);

    /// Update each cell's old_velocity to be the curr_velocity
    for (auto &kv : m_waterGrid) {
        kv.second.oldVelocity = kv.second.currVelocity;
    }
}

Vector3f System::traceParticle(float x, float y, float z, float t) {
    Vector3f vel = getVelocity(Vector3f(x, y, z));
    vel = getVelocity(Vector3f(x + 0.5*t*vel.x(), y + 0.5*t*vel.y(), z + 0.5*t*vel.z()));
    return Vector3f(x, y, z) + t*vel;
}

Vector3i System::getCellIndexFromPoint(Vector3f &pos) {
    return Vector3i{floor(pos.x()), floor(pos.y()), floor(pos.z())};
}

/// Applies the convection term in the Navier-Stokes equation to each cell's velocity
void System::applyConvection(float timeStep) {
    /// Update each cell's curr_velocity in the waterGrid...
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                Cell currCell = m_waterGrid.at(Vector3i{i, j, k});
                Vector3f virtualParticlePos = traceParticle(i + (CELL_DIM/2.f), j + (CELL_DIM/2.f), k + (CELL_DIM/2.f), -timeStep);

                if (!isInBounds(virtualParticlePos.x(), virtualParticlePos.y(), virtualParticlePos.z())) {
                    continue;
                }

                currCell.currVelocity = m_waterGrid.at(getCellIndexFromPoint(virtualParticlePos)).oldVelocity;
            }
        }
    }
}

/// Applies the external force term in the Navier-Stokes equation to each cell's velocity
void System::applyExternalForces(float timeStep) {
    for (auto &kv : m_waterGrid) {
        /// Applies gravity
        kv.second.currVelocity += timeStep * gravity;

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
                m_waterGrid.at(Vector3i(i, j, k)).currVelocity += timeStep * VISCOSITY * Vector3f{u_x, u_y, u_z};
            }
        }
    }
}

void System::initPressureA() {
    int n = WATERGRID_X * WATERGRID_Y * WATERGRID_Z;
    SpMat A(n, n);
    A.setZero();
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                int row_idx = grid2mat(i, j, k);
                std::vector<Vector3i> neighbors = getGridNeighbors(i, j, k);
                for (Vector3i neighbor : neighbors) {
                    A.insert(row_idx, grid2mat(neighbor[0], neighbor[1], neighbor[2])) = 1;
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
    VectorXf b(n, 1);
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                int row_idx = grid2mat(i, j, k);
                float divergence =  getDivergence(i, j, k);
                std::vector<Vector3i> neighbors = getGridNeighbors(i, j, k);
                int ki = 6 - neighbors.size();
                b[row_idx] = ((DENSITY * CELL_DIM) / timeStep) * divergence - ki * ATMOSPHERIC_PRESSURE;
            }
        }
    }
    VectorXf pressure = llt.solve(b);
    return pressure;
}

/// Applies the pressure term in the Navier-Stokes equation to each cell's velocity
void System::applyPressure(float timeStep) {
    VectorXf pressure = calculatePressure(timeStep);
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                int row_idx = grid2mat(i, j, k);
                Vector3f gradient = getGradient(i, j, k, pressure);
                m_waterGrid.at(Vector3i(i, j, k)).currVelocity -= (timeStep / DENSITY * CELL_DIM) * gradient;
            }
        }
    }
}
