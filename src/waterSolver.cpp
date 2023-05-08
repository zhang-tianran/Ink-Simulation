#include "system.h"
#include <math.h>


using namespace Eigen;

double System::updateWaterGrid(double timeToNextRender) {
    double timeStep = calcTimeStep(timeToNextRender);
    updateVelocityField(timeStep);
    return timeStep;
}

double System::calcTimeStep(double timeToNextRender) {
    float maxVelocity = 0;
    for (auto kv : m_waterGrid) {
        if (maxVelocity < kv.second.oldVelocity.norm()) {
            maxVelocity = kv.second.oldVelocity.norm();
        }
    }


//    std::cout << "start" << std::endl;
//    std::cout << maxVelocity << std::endl;

    if (isinf(maxVelocity))
        return MIN_TIMESTEP;

    if (maxVelocity == 0) {
        return MIN_TIMESTEP;
    } else {

        float timeStep = K_CFL * (CELL_DIM / maxVelocity);

        timeStep = std::max(std::min(timeStep, MAX_TIMESTEP), MIN_TIMESTEP);
//        std::cout << timeStep << std::endl;

        assert(timeStep <= MAX_TIMESTEP && timeStep >= MIN_TIMESTEP);
        timeStep = std::min(timeToNextRender, (double)timeStep);
        return timeStep;
    }
}

/// See "Fluid Flow 4 the Rest of Us" Paper for more details
void System::updateVelocityField(float timeStep) {
    /// Error Checking
    #pragma omp parallel for
    for (auto &kv : m_waterGrid) {
        assert(kv.second.oldVelocity == kv.second.currVelocity);
    }

    /// debugging
    checkNanAndInf();

    /// Navier-Stokes equation
    applyConvection(timeStep);
    checkNanAndInf();

    applyExternalForces(timeStep);
    checkNanAndInf();

    applyViscosity(timeStep);
    checkNanAndInf();

    applyPressure(timeStep);
    checkNanAndInf();

    /// Update each cell's old_velocity to be the curr_velocity
    #pragma omp parallel for
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
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                Vector3f virtualParticlePos = traceParticle(i + (CELL_DIM/2.f), j + (CELL_DIM/2.f), k + (CELL_DIM/2.f), -timeStep);

                if (!isInBounds(virtualParticlePos.x(), virtualParticlePos.y(), virtualParticlePos.z())) {
                    continue;
                }

                m_waterGrid[Vector3i(i, j, k)].currVelocity = m_waterGrid.at(getCellIndexFromPoint(virtualParticlePos)).oldVelocity;

                // calculate curl
                m_waterGrid[Vector3i(i, j, k)].curl = getCurl(i, j, k);
            }
        }
    }
}

/// returns a random float [-1, 1]
float zeroOneNoise() {
    float noise = (rand() % 10) / 10.f;
    if (noise > 0.5f) { noise *= -1.f; }
    return noise;
}

Vector3f System::applyWhirlPoolForce(Vector3i index) {
    // Calculate center (independent of y pos)
    Vector3f center(WATERGRID_X/2.F, 0, WATERGRID_Z/2.F);

    // Get unit vector from center axis to index
    Vector3f v(index.x()-center.x(), 0, index.z()-center.z());
    v.normalize();

    // Get perpendicular
    Vector3f whirl = Vector3f(0, 1.f, 0).cross(v);

    // Add some noise
    whirl[0] += zeroOneNoise() * 0.1;
    whirl[1] += zeroOneNoise();
    whirl[2] += zeroOneNoise() * 0.1;

    // Return the cross product
    return whirl * 0.4f;
}

/// Applies the external force term in the Navier-Stokes equation to each cell's velocity
//void System::applyExternalForces(float timeStep) {
//    #pragma omp parallel for collapse(3)
//    for (int i = 0; i < WATERGRID_X; i++) {
//        for (int j = 0; j < WATERGRID_Y; j++) {
//            for (int k = 0; k < WATERGRID_Z; k++) {
//                m_waterGrid[Vector3i(i, j, k)].currVelocity += timeStep * gravity; /// Apply gravity
//                m_waterGrid[Vector3i(i, j, k)].currVelocity += timeStep * applyWhirlPoolForce(Vector3i(i, j, k)); /// Apply whirlpool force
//               // m_waterGrid[Vector3i(i, j, k)].currVelocity += timeStep * getVort(i, j, k); /// Apply vorticity confinement
//            }
//        }
//    }
//}

void System::applyExternalForces(float timeStep) {
<<<<<<< HEAD
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                m_waterGrid[Vector3i(i, j, k)].currVelocity += timeStep * gravity; /// Apply gravity
                m_waterGrid[Vector3i(i, j, k)].currVelocity += timeStep * applyWhirlPoolForce(Vector3i(i, j, k)); /// Apply whirlpool force
                m_waterGrid[Vector3i(i, j, k)].currVelocity += timeStep * getVort(i, j, k); /// Apply vorticity confinement
=======
    std::unordered_set<Vector3i, hash_func> cellsForcesApplied;
    assert(BUFFER_SIZE >=1);
    for (int i = 0; i < m_ink.size(); i++) {
        Vector3i centerCellIndices = getCellIndexFromPoint(m_ink[i].position);
        cellsForcesApplied.insert(centerCellIndices);
        std::vector<Vector3i> neighbors = m_waterGrid[centerCellIndices].neighbors;
        for (int j = 0; j < neighbors.size(); j++) {
            cellsForcesApplied.insert(neighbors[i]);
            if (BUFFER_SIZE >=2) {
                std::vector<Vector3i> neighbors2 = m_waterGrid[neighbors[i]].neighbors;
//                cellsForcesApplied.insert(neighbors2.begin(), neighbors2.end());
                for (int j = 0; j < neighbors2.size(); j++) {
                    cellsForcesApplied.insert(neighbors2[i]);
                    if (BUFFER_SIZE >=3) {
                        std::vector<Vector3i> neighbors3 = m_waterGrid[neighbors2[i]].neighbors;
                        cellsForcesApplied.insert(neighbors3.begin(), neighbors3.end());
                    }
                }
>>>>>>> d482816f9827a9cb8f1b3a17741146b59519f7f7
            }
        }
    }

    for (Vector3i cellIdx: cellsForcesApplied) {
        updateForce(cellIdx, timeStep);
    }
}

void System::updateForce(Vector3i idx, double timeStep){
    m_waterGrid[idx].currVelocity += timeStep * gravity; /// Apply gravity
    m_waterGrid[idx].currVelocity += timeStep * applyWhirlPoolForce(idx); /// Apply whirlpool force
    m_waterGrid[idx].currVelocity += timeStep * getVort(idx[0],idx[1],idx[2]); /// Apply vorticity confinement NOTE THIS WAS SHIFTED IN MAIN
//    auto v = m_waterGrid[idx];
    m_waterGrid[idx].forceApplied = true;
}

Vector3f System::getVort(int i, int j, int k){
    Vector3f curl = m_waterGrid[Vector3i(i, j, k)].curl;
    if (curl == Vector3f(0, 0, 0)) {
        return Vector3f(0, 0, 0);
    }
    Vector3f N = getCurlGradient(i, j, k) / curl.norm();
    Vector3f F_vort = K_VORT * (N.cross(curl));
//    std::cout << "curl" << curl[0] << "," << curl[1] << "," << curl[2] << std::endl;
//    std::cout << "N" << N[0] << "," << N[1] << "," << N[2] << std::endl;
    return F_vort;
}

/// Applies the viscosity term in the Navier-Stokes equation to each cell's velocity
void System::applyViscosity(float timeStep) {
    /// For each cell, apply the viscosity to each cell's currVelocity
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                float u_x = laplacianOperatorOnVelocity(i, j, k, 0);
                float u_y = laplacianOperatorOnVelocity(i, j, k, 1);
                float u_z = laplacianOperatorOnVelocity(i, j, k, 2);
                m_waterGrid[Vector3i(i, j, k)].currVelocity += timeStep * VISCOSITY * Vector3f(u_x, u_y, u_z);
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
                assert(row_idx >= 0 && row_idx < n);
                std::vector<Vector3i> neighbors = getGridNeighbors(i, j, k);
                assert(neighbors.size() <= 6 && neighbors.size() > 2);
                for (Vector3i neighbor : neighbors) {
                    A.insert(row_idx, grid2mat(neighbor[0], neighbor[1], neighbor[2])) = 1;
                }

                A.insert(row_idx, row_idx) = -6;
            }
        }
    }
    llt.compute(A);
    assert(llt.info() == Eigen::Success);
}

/// AP = B (equation 13)
VectorXf System::calculatePressure(float timeStep) {
    int n = WATERGRID_X * WATERGRID_Y * WATERGRID_Z;
    VectorXf b(n, 1);
    #pragma omp parallel for collapse(3)
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
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                Vector3f gradient = getGradient(i, j, k, pressure);
                assert(gradient.norm() < 10000);
                m_waterGrid[Vector3i(i, j, k)].currVelocity -= (timeStep / (DENSITY * CELL_DIM)) * gradient;
            }
        }
    }
}
