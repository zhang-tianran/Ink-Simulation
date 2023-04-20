#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<float> SpMat;

// ============== Global Constants ==============
const int WATERGRID_Z        = 4; /// Water grid width
const int WATERGRID_Y        = 4; /// Water grid height
const int WATERGRID_X        = 4; /// Water grid length
const float CELL_DIM         = 1; /// Cell dimension (is a cube, so length == width == height)

const float DENSITY          = 1; /// Fluid density
const float VISCOSITY        = 1.0016; /// Fluid viscosity. The higher the viscosity, the thicker the liquid.
const float ATMOSPHERIC_PRESSURE = 1; /// Starting number of particles

const int INIT_NUM_PARTICLES = 5; /// Starting number of particles

const Eigen::Vector3f gravity = Eigen::Vector3f(0, -.98, 0);
// ==============================================

typedef struct Cell {
    Eigen::Vector3f oldVelocity;
    Eigen::Vector3f currVelocity;

    // enable printing for debugging
    friend std::ostream& operator<<(std::ostream& strm, const Cell& obj);
} Cell;

typedef struct Particle {
    Eigen::Vector3f position;
    Eigen::Vector3f velocity;
    float opacity;  /// In [0, 1]
    float lifeTime;

    // enable printing for debugging
    friend std::ostream& operator<<(std::ostream& strm, const Particle& obj);
} Particle;

struct hash_func {
    size_t operator()(const Eigen::Vector3i &v) const
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
    const std::vector<Particle>& getInkParticles();

    // enable printing for debugging
    friend std::ostream& operator<<(std::ostream& strm, const System& obj);
private:

    /// Water Grid
    std::unordered_map<Eigen::Vector3i, Cell, hash_func> m_waterGrid;
    void  initWaterGrid();
    float updateWaterGrid();
    float calcTimeStep();
    void  updateVelocityField(float timeStep);
    Eigen::Vector3f traceParticle(float x, float y, float z, float t);
    void  applyConvection(float timeStep);
    void  applyExternalForces(float timeStep);
    float laplacianOperatorOnVelocity(int i, int j, int k, int idx);
    void  applyViscosity(float timeStep);
    Eigen::MatrixXf calculatePressure(float timeStep);
    void  applyPressure(float timeStep);
    int grid2mat(int i, int j, int k) {
        return i * WATERGRID_X + j * WATERGRID_Y + k;
    };
    Eigen::SimplicialLLT<SpMat> llt;
    void initPressureA();

    /// Ink
    std::vector<Particle> m_ink;
    void initParticles();
    void updateParticles(float timeStep);

    /// Getters
    Eigen::Vector3f getGradient(int i, int j, int k, Eigen::VectorXf g);
    float           getDivergence(int i, int j, int k);
    Eigen::Vector3f getVelocity(Eigen::Vector3f pos);
    Eigen::Vector3i getCellIndexFromPoint(Eigen::Vector3f &pos);
    float           getInterpolatedValue(float x, float y, float z, int idx);
    std::vector<Eigen::Vector3i> getGridNeighbors(int i, int j, int k);
    
    /// Boundary Checking: Check if a point (x, y, z) is in bounds of the water grid
    bool isInBounds(float x, float y, float z);
    bool isInBoundsbyIdx(int i, int j, int k);

};


