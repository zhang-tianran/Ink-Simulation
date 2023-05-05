#include "system.h"
#include <stdlib.h>

using namespace Eigen;
using namespace std;

System::System(){}

void System::init() {
    /// Initialize water grid
    initWaterGrid();

    /// Initialize ink particles
    initParticles();
}


/**
 * Initializes a (WATERGRID_X x WATERGRID_Y x WATERGRID_Z) waterGrid
 *               ________________+ (WATERGRID_X+1, WATERGRID_Y+1, WATERGRID_Z+1)
 *             /               / |
 *           /               /   | <---- (WATERGRID_Y)
 *          /--------------/     |
 *          |              |    /
 *          |              |  / <----- (WATERGRID_Z)
 * (0,0,0) +|______________|/
 *           (WATERGRID_X)
 */
void System::initWaterGrid() {
    for (int i = 0; i < WATERGRID_X; i++) {
        for (int j = 0; j < WATERGRID_Y; j++) {
            for (int k = 0; k < WATERGRID_Z; k++) {
                /// Create the cell
                Cell cell {
                    .oldVelocity  = Vector3f(0, 0, 0), // CUSTOMIZABLE
                    .currVelocity = Vector3f(0, 0, 0), // CUSTOMIZABLE
                };

                /// Insert into m_waterGrid
                m_waterGrid.insert({Vector3i{i, j, k}, cell});
            }
        }
    }
    float product = WATERGRID_X*WATERGRID_Z*WATERGRID_Y;
    float size = m_waterGrid.size();
    assert(size == WATERGRID_X*WATERGRID_Z*WATERGRID_Y);
    initPressureA();
}

int nonZeroRand() {
    int ret = rand();
    return (ret != 0) ? ret : 1;
}

/// Returns a random position within the specified ranges
Vector3f getRandPosWithinRange(float minX, float maxX,
                               float minY, float maxY,
                               float minZ, float maxZ) {
    float x = minX + ((maxX - minX) / (nonZeroRand() % 100 + 1));
    assert(!isinf(x) || !isnan(x));
    float y = minY + ((maxY - minY) / (nonZeroRand() % 100 + 1));
    assert(!isinf(y) || !isnan(y));
    float z = minZ + ((maxZ - minZ) / (nonZeroRand() % 100 + 1));
    assert(!isinf(z) || !isnan(z));

    return Vector3f(x, y, z);
}

/// Initializes INIT_NUM_PARTICLES Particle structs
void System::initParticles() {
    m_ink.reserve(INIT_NUM_PARTICLES);
    for (int i = 0; i < INIT_NUM_PARTICLES; i++) {
        /// Create the particle
        Particle particle {
            .position = getRandPosWithinRange(WATERGRID_X/4.f, WATERGRID_X*3/4.f, WATERGRID_Y - 0.1, WATERGRID_Y - 0.1, WATERGRID_Z/4.f, WATERGRID_Z*3/4.f), // CUSTOMIZABLE
            .velocity = Vector3f{0, 0, 0}, // CUSTOMIZABLE
            .opacity  = 1.f,
            .lifeTime = 5.f // CUSTOMIZABLE
        };

        /// Insert into m_ink
        m_ink.push_back(particle);
    }
    assert(m_ink.size() == INIT_NUM_PARTICLES);
}

/************************** GETTERS ************************************/
const std::vector<Particle>& System::getInkParticles() {
    return m_ink;
}
