#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat;

// ============== Global Constants ==============
const int WATERGRID_X        = 8; /// Water grid length
const int WATERGRID_Y        = 8; /// Water grid height
const int WATERGRID_Z        = 8; /// Water grid width
const double CELL_DIM         = 1; /// Cell dimension (is a cube, so length == width == height)

const double DENSITY          = 1; /// Fluid density

const double K_VORT           = 1; /// strength of vorticity

const double VISCOSITY        = 1.0016; /// 1.0016  /// Fluid viscosity. The higher the viscosity, the thicker the liquid.
const double ATMOSPHERIC_PRESSURE = 1; /// Starting number of particles

const int INIT_NUM_PARTICLES = 5000; /// Starting number of particles

const Eigen::Vector3d gravity = Eigen::Vector3d(0, -0.58, 0);
//const Eigen::Vector3d gravity = Eigen::Vector3d(0, -0.98, 0);

const double K_CFL = 0.2f;
const double MIN_TIMESTEP = 0.01f;
const double MAX_TIMESTEP = 1.f;
// ==============================================

enum CellBFECCField {
    OLDVELOCITY,
    USTARFORWARD,
    USTAR,
    USQUIGGLY,
    CURRVELOCITY
};

typedef struct Cell {
    Eigen::Vector3d oldVelocity;
    Eigen::Vector3d currVelocity;

    /// For BFECC
    Eigen::Vector3d uStarForward;
    Eigen::Vector3d uStar;
    Eigen::Vector3d uSquiggly;

    Eigen::Vector3d curl;

    bool forceApplied;
    std::vector<Eigen::Vector3i> neighbors;

    // enable printing for debugging
    friend std::ostream& operator<<(std::ostream& strm, const Cell& obj);
} Cell;

typedef struct Particle {
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    double opacity;  /// In [0, 1]
    double lifeTime;

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

    // check if particles or watergrid values for pos/vel have inf or nans
    void checkNanAndInf();

    // enable printing for debugging
    friend std::ostream& operator<<(std::ostream& strm, const System& obj);
private:

    /// Water Grid
    std::unordered_map<Eigen::Vector3i, Cell, hash_func> m_waterGrid;
    void  initWaterGrid();
    double updateWaterGrid();
    double calcTimeStep();
    void  updateVelocityField(double timeStep);
    Eigen::Vector3d traceParticle(double x, double y, double z, double t);
     Eigen::Vector3d traceParticle(double x, double y, double z, double t, CellBFECCField field);
    void  applyConvection(double timeStep, CellBFECCField field);
    void  applyExternalForces(double timeStep);
    void  updateForce(Eigen::Vector3i idx, double timeStep);
    Eigen::Vector3d getVort(Eigen::Vector3i idx);
    void  applyViscosity(double timeStep);
    Eigen::VectorXd calculatePressure(double timeStep);
    void  applyPressure(double timeStep);
    Eigen::SparseLU<SpMat> llt;
    void initPressureA();

    // vorticity
    void applyVorticity(double timestep);

    void applyBFECC(double timeStep);

    Eigen::Vector3d applyWhirlPoolForce(Eigen::Vector3i index);

    int grid2mat(int i, int j, int k) {
        return (i * WATERGRID_Z * WATERGRID_Y) + (j * WATERGRID_X) + k;
    };

    /// Ink
    std::vector<Particle> m_ink;
    void initParticles();
    void updateParticles(double timeStep);

    /// Getters
    Eigen::Vector3d getGradient(int i, int j, int k, Eigen::VectorXd g);
    double           getDivergence(int i, int j, int k);
    Eigen::Vector3d getCurl(int i, int j, int k);
    Eigen::Vector3d getCurlGradient(int i, int j, int k);
    double           laplacianOperatorOnVelocity(int i, int j, int k, int idx);
    Eigen::Vector3d getVelocity(Eigen::Vector3d pos);
    Eigen::Vector3d getVelocity(Eigen::Vector3d pos, CellBFECCField field);
    Eigen::Vector3i getCellIndexFromPoint(Eigen::Vector3d &pos);
    double           getInterpolatedValue(double x, double y, double z, int idx);
    double           getInterpolatedValue(double x, double y, double z, int idx, CellBFECCField field);
    std::vector<Eigen::Vector3i> getGridNeighbors(int i, int j, int k);
    Eigen::Vector3d getVelocityFromField(Eigen::Vector3i pos, CellBFECCField field);
    
    /// Boundary Checking: Check if a point (x, y, z) is in bounds of the water grid
    bool isInBounds(double x, double y, double z);
    bool isInBoundsbyIdx(int i, int j, int k);

    /// DEBUGGING
   bool hasNan(Eigen::Vector3d v);
   bool hasInf(Eigen::Vector3d v);
};


