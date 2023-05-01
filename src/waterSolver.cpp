#include "system.h"
#include <math.h>

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

    if (isinf(maxVelocity))
        return MIN_TIMESTEP;

    if (maxVelocity == 0) {
        return MIN_TIMESTEP;
    } else {

        float timeStep = K_CFL * (CELL_DIM / maxVelocity);

        timeStep = std::max(std::min(timeStep, MAX_TIMESTEP), MIN_TIMESTEP);

        assert(timeStep <= MAX_TIMESTEP && timeStep >= MIN_TIMESTEP);
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
#ifndef QT_NO_DEBUG
    checkNanAndInf();
#endif

    /// Navier-Stokes equation
    applyBFECC(timeStep);
#ifndef QT_NO_DEBUG
    checkNanAndInf();
#endif

    applyExternalForces(timeStep);
#ifndef QT_NO_DEBUG
    checkNanAndInf();
#endif

    applyViscosity(timeStep);
#ifndef QT_NO_DEBUG
    checkNanAndInf();
#endif

    applyPressure(timeStep);
#ifndef QT_NO_DEBUG
    checkNanAndInf();
#endif

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

Eigen::Vector3f System::traceParticle(float x, float y, float z, float t, CellBFECCField field) {
    Vector3f vel = getVelocity(Vector3f(x, y, z), field);
    vel = getVelocity(Vector3f(x + 0.5*t*vel.x(), y + 0.5*t*vel.y(), z + 0.5*t*vel.z()), field);
    return Vector3f(x, y, z) + t*vel;
}

/// Applies the convection term in the Navier-Stokes equation to each cell's velocity
void System::applyConvection(float timeStep, CellBFECCField field) {
    /// Update each cell's curr_velocity in the waterGrid...
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                /// Trace particle, fields are input b/c we want to get the previous iterations velocity (i.e. vel when we want to find USQUIGGLY will be USTAR)
                Vector3f virtualParticlePos; //= traceParticle(i + (CELL_DIM/2.f), j + (CELL_DIM/2.f), k + (CELL_DIM/2.f), -timeStep, field);
                switch(field) {
                    case OLDVELOCITY:
                        virtualParticlePos = traceParticle(i + (CELL_DIM/2.f), j + (CELL_DIM/2.f), k + (CELL_DIM/2.f), -timeStep, field);
                        break;
                    case USTARFORWARD:
                        virtualParticlePos = traceParticle(i + (CELL_DIM/2.f), j + (CELL_DIM/2.f), k + (CELL_DIM/2.f), -timeStep, CURRVELOCITY);
                        break;
                    case USTAR:
                        virtualParticlePos = traceParticle(i + (CELL_DIM/2.f), j + (CELL_DIM/2.f), k + (CELL_DIM/2.f), -timeStep, USTARFORWARD);
                        break;
                    case CURRVELOCITY:
                        virtualParticlePos = traceParticle(i + (CELL_DIM/2.f), j + (CELL_DIM/2.f), k + (CELL_DIM/2.f), -timeStep, USQUIGGLY);
                        break;
                    case USQUIGGLY:
                        virtualParticlePos = traceParticle(i + (CELL_DIM/2.f), j + (CELL_DIM/2.f), k + (CELL_DIM/2.f), -timeStep, USTAR);
                        break;
                }

                if (!isInBounds(virtualParticlePos.x(), virtualParticlePos.y(), virtualParticlePos.z())) {
                    continue;
                }

                /// Updating fields
                Cell *cell = &m_waterGrid.at(getCellIndexFromPoint(virtualParticlePos));
                switch(field) {
                    case OLDVELOCITY:
                        m_waterGrid[Vector3i(i, j, k)].currVelocity = cell->oldVelocity;
                        break;
                    case USTARFORWARD:
                        m_waterGrid[Vector3i(i, j, k)].uStarForward = cell->currVelocity;
                        break;
                    case USTAR:
                        m_waterGrid[Vector3i(i, j, k)].uStar        = cell->uStarForward;
                        break;
                    case CURRVELOCITY:
                        m_waterGrid[Vector3i(i, j, k)].currVelocity = cell->uSquiggly;
                        break;
                    case USQUIGGLY:
                        Vector3f error = (cell->uStar - cell->oldVelocity) / 2.f;
                        m_waterGrid[Vector3i(i, j, k)].uSquiggly    = cell->oldVelocity - error;
                        break;
                }

                /// Calculate curl (just once)
                if (field == CURRVELOCITY) {
                    m_waterGrid[Vector3i(i, j, k)].curl = getCurl(i, j, k);
                }
            }
        }
    }
}

/// Applies BFECC algorithm to convection term in Navier-Stokes equation
void System::applyBFECC(float timeStep) {
    /// 1) apply backwards particle trace to get u*_{n+1}
    applyConvection(-timeStep, USTARFORWARD); /// forward in time update
    
    /// 2) find reverse of backwards particle trace from u*_{n+1} to get u*_{n}
    applyConvection(timeStep, USTAR); /// backward in time update
    
    /// 3) starting from part 3 - do a backwards particle trace for u-squiggle
    applyConvection(-timeStep, USQUIGGLY); /// forward in time update

    /// 4) starting from part 4 - do a backwards particle trace for currVelocity
    applyConvection(-timeStep, CURRVELOCITY); /// forward in time update
}

/// returns a random float [-1, 1]
float zeroOneNoise() {
    float noise = (rand() % 10) / 10.f;
    if (noise > 0.5f) { noise *= -1.f; }
    return noise;
}

/// Calculates an whirlpool external force
Vector3f System::applyWhirlPoolForce(Vector3i index) {
    /// Calculate center (independent of y pos)
    Vector3f center(WATERGRID_X/2.F, 0, WATERGRID_Z/2.F);

    /// Get unit vector from center axis to index
    Vector3f v(index.x()-center.x(), 0, index.z()-center.z());
    v.normalize();

    /// Get perpendicular
    Vector3f whirl = Vector3f(0, 1.f, 0).cross(v);

    /// Add some noise
    whirl[0] += zeroOneNoise() * 0.1;
    whirl[1] += zeroOneNoise();
    whirl[2] += zeroOneNoise() * 0.1;

    /// Return the cross product
    return whirl * 0.4f;
}

/// Applies the external force term in the Navier-Stokes equation to each cell's velocity
void System::applyExternalForces(float timeStep) {
    std::vector<Vector3i> cellsForcesApplied;
    for (int i = 0; i < m_ink.size(); i++) {
        Particle currParticle = m_ink[i];
        Vector3i centerCellIndices = getCellIndexFromPoint(currParticle.position);
        Cell centerCell = m_waterGrid[centerCellIndices];
        if (centerCell.forceApplied == false) {
            updateForce(centerCellIndices, timeStep);
            cellsForcesApplied.push_back(centerCellIndices);
        }
        std::vector<Vector3i> neighbors = centerCell.neighbors;
        for (int j = 0; j < neighbors.size(); j++) {
            if (m_waterGrid[neighbors[j]].forceApplied == false) {
                updateForce(neighbors[j], timeStep);
                cellsForcesApplied.push_back(neighbors[j]);
            }
        }
    }

    #pragma omp parallel for
    for (Vector3i cellIdx: cellsForcesApplied) {
        m_waterGrid[cellIdx].forceApplied = false;
    }
}

void System::updateForce(Vector3i idx, float timeStep){

    m_waterGrid[idx].currVelocity += timeStep * gravity; /// Apply gravity
//    m_waterGrid[idx].currVelocity += timeStep * applyWhirlPoolForce(idx); /// Apply whirlpool force
    m_waterGrid[idx].currVelocity += timeStep * getVort(idx); /// Apply vorticity confinement
    auto v = m_waterGrid[idx];
    m_waterGrid[idx].forceApplied = true;
}

Vector3f System::getVort(Vector3i idx){
    Vector3f curl = m_waterGrid[idx].curl;
    if (curl == Vector3f(0, 0, 0)) {
        return Vector3f(0, 0, 0);
    }
    Vector3f N = getCurlGradient(idx[0], idx[1], idx[2]) / curl.norm();
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
                std::vector<Vector3i> neighbors = m_waterGrid[Vector3i(i, j, k)].neighbors;
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
                int ki = 6 - m_waterGrid[Vector3i(i, j, k)].neighbors.size();
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
