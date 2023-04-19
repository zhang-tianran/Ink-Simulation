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
 * Initializes a (WATERGRID_L x WATERGRID_W x WATERGRID_H) waterGrid
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
    assert(m_waterGrid.size() == WATERGRID_L*WATERGRID_W*WATERGRID_H);
}

/// Returns a random position within the specified ranges
Vector3f getRandPosWithinRange(float minX, float maxX,
                               float minY, float maxY,
                               float minZ, float maxZ) {
    return Vector3f{
        minX + ((maxX - minX) / (rand() % 1000)),
        minY + ((maxY - minY) / (rand() % 1000)),
        minZ + ((maxZ - minZ) / (rand() % 1000)),
    };
}

/// Initializes INIT_NUM_PARTICLES Particle structs
void System::initParticles() {
    m_ink.reserve(INIT_NUM_PARTICLES);
    for (int i = 0; i < INIT_NUM_PARTICLES; i++) {
        /// Create the particle
        Particle particle {
            .position = getRandPosWithinRange(WATERGRID_L/4.f, WATERGRID_L*3/4.f, WATERGRID_H - 0.001, WATERGRID_H - 0.001, WATERGRID_W/4.f, WATERGRID_W*3/4.f), // CUSTOMIZABLE
            .velocity = Vector3f{0, -5, 0}, // CUSTOMIZABLE
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

/************************** PRINTING UTILS *****************************/

/// print a cell
ostream& operator<<(ostream& strm, const Cell& obj) {
    strm << "\tvelocity: (" << obj.velocity.x() << ", ";
    strm << obj.velocity.y() << ", " << obj.velocity.z() << ")\n";
    strm << "\tpressure: " << obj.pressure;
    return strm;
}

/// print a particle
ostream& operator<<(ostream& strm, const Particle& obj) {
    strm << "Particle: \n";
    strm << "\tpos: (" << obj.position.x() << ", ";
    strm << obj.position.y() << ", " << obj.position.z() << ")\n";
    strm << "\tvelocity: (" << obj.velocity.x() << ", ";
    strm << obj.velocity.y() << ", " << obj.velocity.z() << ")\n";
    strm << "\topacity: " << obj.opacity;
    strm << "\tlifetime: " << obj.lifeTime;
    return strm;
}

/// print the whole system
ostream& operator<<(ostream& strm, const System& obj) {
    strm << "********* PRINTING SYSTEM ***********\n";
    strm << "********* PRINTING CELLS ***********\n";
    for (auto& [k, v] : obj.m_waterGrid) {
        strm << "Cell: \n";
        strm << "\tpos in hashmap: " << k.x() << ", ";
        strm << k.x() << ", " << k.z() << ")\n";
        strm << v << endl;
    }

    strm << "********* PRINTING PARTICLES ***********\n";
    for (auto& el : obj.m_ink) {
        strm << el << endl;
    }
    return strm;
}
