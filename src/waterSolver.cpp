#include "system.h"
#include <math.h>
#include <unordered_set>

using namespace Eigen;

#define DO_BFECC false

double System::updateWaterGrid() {
    double timeStep = calcTimeStep();
    updateVelocityField(timeStep);
    return timeStep;
}


double System::calcTimeStep() {
    double maxVelocity = 0;
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

        double timeStep = K_CFL * (CELL_DIM / maxVelocity);

        timeStep = std::max(std::min(timeStep, MAX_TIMESTEP), MIN_TIMESTEP);

        assert(timeStep <= MAX_TIMESTEP && timeStep >= MIN_TIMESTEP);
        return timeStep;
    }
}

/// See "Fluid Flow 4 the Rest of Us" Paper for more details
void System::updateVelocityField(double timeStep) {
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

    /// Update each cell's old_velocity to be the curr_velocity
    #pragma omp parallel for
    for (auto &kv : m_waterGrid) {
        kv.second.oldVelocity = kv.second.currVelocity;
    }

    applyViscosity(timeStep);
#ifndef QT_NO_DEBUG
    checkNanAndInf();
#endif

    /// Update each cell's old_velocity to be the curr_velocity
    #pragma omp parallel for
    for (auto &kv : m_waterGrid) {
        kv.second.oldVelocity = kv.second.currVelocity;
    }

   applyPressure(timeStep);
#ifndef QT_NO_DEBUG
    checkNanAndInf();
#endif

    /// Update each cell's old_velocity to be the curr_velocity
    #pragma omp parallel for
    for (auto &kv : m_waterGrid) {
        kv.second.oldVelocity = kv.second.currVelocity;
    }

//    applyVorticity(timeStep);
    checkNanAndInf();

    /// Update each cell's old_velocity to be the curr_velocity
    #pragma omp parallel for
    for (auto &kv : m_waterGrid) {
        kv.second.oldVelocity = kv.second.currVelocity;
    }
}

Vector3d System::traceParticle(double x, double y, double z, double t) {
    Vector3d vel = getVelocity(Vector3d(x, y, z));
    vel = getVelocity(Vector3d(x + 0.5*t*vel.x(), y + 0.5*t*vel.y(), z + 0.5*t*vel.z()));
    return Vector3d(x, y, z) + t*vel;
}

Eigen::Vector3d System::traceParticle(double x, double y, double z, double t, CellBFECCField field) {
    Vector3d vel = getVelocity(Vector3d(x, y, z), field);
    vel = getVelocity(Vector3d(x + 0.5*t*vel.x(), y + 0.5*t*vel.y(), z + 0.5*t*vel.z()), field);
    return Vector3d(x, y, z) + t*vel;
}

/// Applies the convection term in the Navier-Stokes equation to each cell's velocity
void System::applyConvection(double timeStep, CellBFECCField field) {
    /// Update each cell's curr_velocity in the waterGrid...
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                /// Trace particle, fields are input b/c we want to get the previous iterations velocity (i.e. vel when we want to find USQUIGGLY will be USTAR)
                Vector3d virtualParticlePosX;
                Vector3d virtualParticlePosY;
                Vector3d virtualParticlePosZ;
                Vector3d virtualParticleVelocity{0, 0, 0};
                switch(field) {
                    case OLDVELOCITY: /// Note: This case never gets called/never happens
                        break;
                    case USTARFORWARD: /// 1) Starting from oldVelocity, calculate uStarForward
                        virtualParticlePosX = traceParticle(i +  CELL_DIM,      j + (CELL_DIM/2.f), k + (CELL_DIM/2.f), timeStep, OLDVELOCITY);
                        virtualParticlePosY = traceParticle(i + (CELL_DIM/2.f), j +  CELL_DIM,      k + (CELL_DIM/2.f), timeStep, OLDVELOCITY);
                        virtualParticlePosZ = traceParticle(i + (CELL_DIM/2.f), j + (CELL_DIM/2.f), k +  CELL_DIM,      timeStep, OLDVELOCITY);

                        if (isInBounds(virtualParticlePosX.x(), virtualParticlePosX.y(), virtualParticlePosX.z())) {
                            virtualParticleVelocity.x() = m_waterGrid[getCellIndexFromPoint(virtualParticlePosX)].oldVelocity.x();
                        }
                        if (isInBounds(virtualParticlePosY.x(), virtualParticlePosY.y(), virtualParticlePosY.z())) {
                            virtualParticleVelocity.y() = m_waterGrid[getCellIndexFromPoint(virtualParticlePosY)].oldVelocity.y();
                        }
                        if (isInBounds(virtualParticlePosZ.x(), virtualParticlePosZ.y(), virtualParticlePosZ.z())) {
                            virtualParticleVelocity.z() = m_waterGrid[getCellIndexFromPoint(virtualParticlePosZ)].oldVelocity.z();
                        }
                        break;
                    case USTAR: /// 2) Starting from uStarForward, calculate uStar
                        virtualParticlePosX = traceParticle(i +  CELL_DIM,      j + (CELL_DIM/2.f), k + (CELL_DIM/2.f), timeStep, USTARFORWARD);
                        virtualParticlePosY = traceParticle(i + (CELL_DIM/2.f), j +  CELL_DIM,      k + (CELL_DIM/2.f), timeStep, USTARFORWARD);
                        virtualParticlePosZ = traceParticle(i + (CELL_DIM/2.f), j + (CELL_DIM/2.f), k +  CELL_DIM,      timeStep, USTARFORWARD);

                        if (isInBounds(virtualParticlePosX.x(), virtualParticlePosX.y(), virtualParticlePosX.z())) {
                            virtualParticleVelocity.x() = m_waterGrid[getCellIndexFromPoint(virtualParticlePosX)].uStarForward.x();
                        }
                        if (isInBounds(virtualParticlePosY.x(), virtualParticlePosY.y(), virtualParticlePosY.z())) {
                            virtualParticleVelocity.y() = m_waterGrid[getCellIndexFromPoint(virtualParticlePosY)].uStarForward.y();
                        }
                        if (isInBounds(virtualParticlePosZ.x(), virtualParticlePosZ.y(), virtualParticlePosZ.z())) {
                            virtualParticleVelocity.z() = m_waterGrid[getCellIndexFromPoint(virtualParticlePosZ)].uStarForward.z();
                        }
                        break;
                    case CURRVELOCITY: /// 4) Starting from uSquiggly, calculate currVelocity
                        virtualParticlePosX = traceParticle(i +  CELL_DIM,      j + (CELL_DIM/2.f), k + (CELL_DIM/2.f), timeStep, USQUIGGLY);
                        virtualParticlePosY = traceParticle(i + (CELL_DIM/2.f), j +  CELL_DIM,      k + (CELL_DIM/2.f), timeStep, USQUIGGLY);
                        virtualParticlePosZ = traceParticle(i + (CELL_DIM/2.f), j + (CELL_DIM/2.f), k +  CELL_DIM,      timeStep, USQUIGGLY);

                        if (isInBounds(virtualParticlePosX.x(), virtualParticlePosX.y(), virtualParticlePosX.z())) {
                            virtualParticleVelocity.x() = m_waterGrid[getCellIndexFromPoint(virtualParticlePosX)].uSquiggly.x();
                        }
                        if (isInBounds(virtualParticlePosY.x(), virtualParticlePosY.y(), virtualParticlePosY.z())) {
                            virtualParticleVelocity.y() = m_waterGrid[getCellIndexFromPoint(virtualParticlePosY)].uSquiggly.y();
                        }
                        if (isInBounds(virtualParticlePosZ.x(), virtualParticlePosZ.y(), virtualParticlePosZ.z())) {
                            virtualParticleVelocity.z() = m_waterGrid[getCellIndexFromPoint(virtualParticlePosZ)].uSquiggly.z();
                        }
                        break;
                    case USQUIGGLY: /// 3) Calculate error on the cell we are currently on (not using particle trace on this one aka not actually using virtualParticlePos here)
                        break;
                }

                /// Updating fields
                switch(field) {
                    case OLDVELOCITY: /// Note: This case never gets called/never happens
                        break;
                    case USTARFORWARD: /// 1) Starting from oldVelocity, calculate uStarForward
                        if (DO_BFECC)
                            m_waterGrid[Vector3i(i, j, k)].uStarForward = virtualParticleVelocity;
                        else
                            m_waterGrid[Vector3i(i, j, k)].currVelocity = virtualParticleVelocity;
                        break;
                    case USTAR: /// 2) Starting from uStarForward, calculate uStar
                        m_waterGrid[Vector3i(i, j, k)].uStar        = virtualParticleVelocity;
                        break;
                    case CURRVELOCITY: /// 4) Starting uStar, calculate currVelocity
                        m_waterGrid[Vector3i(i, j, k)].currVelocity = virtualParticleVelocity;
                        break;
                    case USQUIGGLY: /// 3) Calculate error on the cell we are currently on (not using particle trace on this one)
                        Vector3d error = (m_waterGrid[Vector3i(i, j, k)].uStar - m_waterGrid[Vector3i(i, j, k)].oldVelocity) / 2.f;
                        m_waterGrid[Vector3i(i, j, k)].uSquiggly    = m_waterGrid[Vector3i(i, j, k)].oldVelocity - error;
                        break;
                }
            }
        }
    }
}

/// Applies BFECC algorithm to convection term in Navier-Stokes equation
void System::applyBFECC(double timeStep) {
    /// 1) apply backwards particle trace to get u*_{n+1}
    applyConvection(-timeStep, USTARFORWARD); /// forward in time update
    
    if (DO_BFECC) {
        /// 2) find reverse of backwards particle trace from u*_{n+1} to get u*_{n}
        applyConvection(timeStep, USTAR); /// backward in time update

        /// 3) starting from part 3 - estimate error for usquiggle
        applyConvection(-timeStep, USQUIGGLY); /// We're not actually moving in any direction for time step

        /// 4) starting from part 4 - do a backwards particle trace for currVelocity
        applyConvection(-timeStep, CURRVELOCITY); /// forward in time update
    }

}

/// returns a random double [-1, 1]
double zeroOneNoise() {
    double noise = (rand() % 10) / 10.f;
    if (noise > 0.5f) { noise *= -1.f; }
    return noise;
}

/// Calculates an whirlpool external force
Vector3d System::applyWhirlPoolForce(Vector3i index) {
    /// Calculate center (independent of y pos)
    Vector3d center(WATERGRID_X/2.F, 0, WATERGRID_Z/2.F);

    /// Get unit vector from center axis to index
    Vector3d v(index.x()-center.x(), 0, index.z()-center.z());
    v.normalize();

    /// Get perpendicular
    Vector3d whirl = Vector3d(0, 1.f, 0).cross(v);

    /// Add some noise
    whirl[0] += zeroOneNoise() * 0.1;
    whirl[1] += zeroOneNoise();
    whirl[2] += zeroOneNoise() * 0.1;

    /// Return the cross product
    return whirl * 0.4f;
}

///// Applies the external force term in the Navier-Stokes equation to each cell's velocity
//void System::applyExternalForces(float timeStep) {
//    std::vector<Vector3i> cellsForcesApplied;
//    for (int i = 0; i < m_ink.size(); i++) {
//        Particle currParticle = m_ink[i];
//        Vector3i centerCellIndices = getCellIndexFromPoint(currParticle.position);
//        Cell centerCell = m_waterGrid[centerCellIndices];
//        if (centerCell.forceApplied == false) {
//            updateForce(centerCellIndices, timeStep);
//            cellsForcesApplied.push_back(centerCellIndices);
//        }
//        std::vector<Vector3i> neighbors = centerCell.neighbors;
//        for (int j = 0; j < neighbors.size(); j++) {
//            if (m_waterGrid[neighbors[j]].forceApplied == false) {
//                updateForce(neighbors[j], timeStep);
//                cellsForcesApplied.push_back(neighbors[j]);
//            }
//        }
//    }

//    #pragma omp parallel for
//    for (Vector3i cellIdx: cellsForcesApplied) {
//        m_waterGrid[cellIdx].forceApplied = false;
//    }
//}

/// Applies the external force term in the Navier-Stokes equation to each cell's velocity
void System::applyExternalForces(double timeStep) {
    std::unordered_set<Vector3i, hash_func> cellsForcesApplied;
    for (int i = 0; i < m_ink.size(); i++) {
        Vector3i centerCellIndices = getCellIndexFromPoint(m_ink[i].position);
        cellsForcesApplied.insert(centerCellIndices);
        std::vector<Vector3i> neighbors = m_waterGrid[centerCellIndices].neighbors;
        for (int j = 0; j < neighbors.size(); j++) {
            cellsForcesApplied.insert(neighbors[j]);
            std::vector<Vector3i> neighbors2 = m_waterGrid[neighbors[j]].neighbors;
            cellsForcesApplied.insert(neighbors2.begin(), neighbors2.end());
        }
    }
    for (Vector3i cellIdx: cellsForcesApplied) {
        updateForce(cellIdx, timeStep);
    }
}

void System::updateForce(Vector3i idx, double timeStep){

    m_waterGrid[idx].currVelocity += timeStep * gravity; /// Apply gravity
//    m_waterGrid[idx].currVelocity += timeStep * applyWhirlPoolForce(idx); /// Apply whirlpool force
    auto v = m_waterGrid[idx];
    m_waterGrid[idx].forceApplied = true;
}

/// Applies the viscosity term in the Navier-Stokes equation to each cell's velocity
void System::applyViscosity(double timeStep) {
    /// For each cell, apply the viscosity to each cell's currVelocity
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                double u_x = laplacianOperatorOnVelocity(i, j, k, 0);
                double u_y = laplacianOperatorOnVelocity(i, j, k, 1);
                double u_z = laplacianOperatorOnVelocity(i, j, k, 2);
                m_waterGrid[Vector3i(i, j, k)].currVelocity += timeStep * VISCOSITY * Vector3d(u_x, u_y, u_z);
            }
        }
    }
}


void System::applyVorticity(double timestep) {
    /// For each cell, calculate the curl
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                m_waterGrid[Vector3i(i, j, k)].curl = getCurl(i, j, k);
            }
        }
    }

    /// For each cell, apply the vorticity confinement force
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                m_waterGrid[Vector3i(i, j, k)].currVelocity += timestep * getVort(Vector3i(i, j, k));
            }
        }
    }
}


Vector3d System::getVort(Vector3i idx){
    Vector3d curl = m_waterGrid[idx].curl;
    if (curl == Vector3d(0, 0, 0)) {
        return Vector3d(0, 0, 0);
    }
    Vector3d N = getCurlGradient(idx[0], idx[1], idx[2]) / curl.norm();
    Vector3d F_vort = K_VORT * (N.cross(curl));
    return F_vort;
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
VectorXd System::calculatePressure(double timeStep) {
    int n = WATERGRID_X * WATERGRID_Y * WATERGRID_Z;
    VectorXd b(n, 1);
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                int row_idx = grid2mat(i, j, k);
                double divergence =  getDivergence(i, j, k);
                int ki = 6 - m_waterGrid[Vector3i(i, j, k)].neighbors.size();
                b[row_idx] = ((DENSITY * CELL_DIM) / timeStep) * divergence - ki * ATMOSPHERIC_PRESSURE;
            }
        }
    }
    VectorXd pressure = llt.solve(b);
    return pressure;
}

/// Applies the pressure term in the Navier-Stokes equation to each cell's velocity
void System::applyPressure(double timeStep) {
    VectorXd pressure = calculatePressure(timeStep);
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                Vector3d gradient = getGradient(i, j, k, pressure);
                assert(gradient.norm() < 10000);
                m_waterGrid[Vector3i(i, j, k)].currVelocity -= (timeStep / (DENSITY * CELL_DIM)) * gradient;
            }
        }
    }
}
