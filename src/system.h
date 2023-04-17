#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>

// ============== Global Constants ==============
const int WATERGRID_W        = 4; /// Water grid width
const int WATERGRID_H        = 4; /// Water grid height
const int WATERGRID_L        = 4; /// Water grid length
const float CELL_DIM         = 1; /// Cell dimension (is a cube, so length == width == height)

const float DENSITY          = 1; /// Fluid density
const float VISCOSITY        = 1.0016; /// Fluid viscosity. The higher the viscosity, the thicker the liquid.

const int INIT_NUM_PARTICLES = 5; /// Starting number of particles
// ==============================================

typedef struct Cell {
    Eigen::Vector3f velocity;
    float pressure;
} Cell;

typedef struct Particle {
    Eigen::Vector3f position;
    Eigen::Vector3f velocity;
    float opacity;  /// In [0, 1]
    float lifeTime;
} Particle;

struct hash_func {
    size_t operator()(const Eigen::Vector3f &v) const
    {
        assert(v.x()>=0 && v.y()>=0 && v.z()>=0);
        return 541 * v.x() + 79 * v.y() + 31 * v.z();
    }
};

class System
{
public:
    System();

    void init();
    void solve();

private:
    /// Water Grid
    std::unordered_map<Eigen::Vector3f, Cell, hash_func> m_waterGrid;
    void  initWaterGrid();
    float updateWaterGrid();
    float calcTimeStep();
    void  updateVelocityField(float timeStep);
    void  applyConvection(float timeStep);
    void  applyExternalForces(float timeStep);
    void  applyViscosity(float timeStep);
    void  calculatePressure();
    void  applyPressure(float timeStep);

    /// Ink
    std::vector<Particle> m_ink;
    void initParticles();
    void updateParticles(float timeStep);
};


