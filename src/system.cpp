#include "system.h"

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
 * Initializes a (WATERGRID_W x WATERGRID_H x WATERGRID_L) waterGrid
 *               ________________+ (WATERGRID_L, WATERGRID_W, WATERGRID_H)
 *             /               / |
 *           /               /   | <---- (WATERGRID_H)
 *          /--------------/     |
 *          |              |    /
 *          |              |  / <----- (WATERGRID_W)
 * (0,0,0) +|______________|/
 *           (WATERGRID_L)
 */
void System::initWaterGrid() {
    for (int l = 0; l < WATERGRID_L; l++) {
        for (int w = 0; w < WATERGRID_W; w++) {
            for (int h = 0; h < WATERGRID_H; h++) {
                /// Create the cell
                Cell cell {
                    .velocity = Vector3f{5, 0, 0}, // CUSTOMIZABLE
                    .pressure = 0
                };

                /// Insert into m_waterGrid
                m_waterGrid.insert({Vector3f{l, w, h}, cell});
            }
        }
    }
}

/// Returns a random position within the specified ranges
Vector3f getRandPosWithinRange(float minX, float maxX,
                               float minY, float maxY,
                               float minZ, float maxZ) {
    return Vector3f{
        minX + ((maxX - minX) / (arc4random() % 1000)),
        minY + ((maxY - minY) / (arc4random() % 1000)),
        minZ + ((maxZ - minZ) / (arc4random() % 1000)),
    };
}

/// Initializes INIT_NUM_PARTICLES Particle structs
void System::initParticles() {
    m_ink.reserve(INIT_NUM_PARTICLES);
    for (int i = 0; i < INIT_NUM_PARTICLES; i++) {
        /// Create the particle
        Particle particle {
            .position = getRandPosWithinRange(WATERGRID_L/4.f, WATERGRID_L*3/4.f, 0.f, 0.f, WATERGRID_W/4.f, WATERGRID_W*3/4.f), // CUSTOMIZABLE
            .velocity = Vector3f{0, -5, 0}, // CUSTOMIZABLE
            .opacity  = 1.f,
            .lifeTime = 5.f
        };

        /// Insert into m_ink
        m_ink.push_back(particle);
    }
}


