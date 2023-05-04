#include "system.h"
using namespace Eigen;
using namespace std;

#define EPSILON .005
Vector3i System::getCellIndexFromPoint(Vector3d &pos) {
    return Vector3i{floor(pos.x() / CELL_DIM), floor(pos.y() / CELL_DIM), floor(pos.z() / CELL_DIM)};
}

Vector3d System::getGradient(int i, int j, int k, VectorXd g){
    int row_idx = grid2mat(i, j, k);
    double val = g[row_idx];
    Vector3d gradient(val, val, val);
    if (isInBoundsbyIdx(i - 1, j, k)) {
        gradient[0] -= g[grid2mat(i - 1, j, k)];
    } else {
        gradient[0] -= 1;
    }

    if (isInBoundsbyIdx(i, j - 1, k)) {
        gradient[1] -= g[grid2mat(i, j - 1, k)];
    } else {
        gradient[1] -= 1;
    }

    if (isInBoundsbyIdx(i, j, k - 1)) {
        gradient[2] -= g[grid2mat(i, j, k - 1)];
    } else {
        gradient[2] -= 1;
    }

    if (gradient.norm() > 1000) {
//        std::cout<<"GRADIENT" << gradient<<std::endl;
//        std::cout<<"PRESSURE" << g<<std::endl;
        assert(gradient.norm() < 1000);
    }
    return gradient;
}

double System::getDivergence(int i, int j, int k){
    Vector3d vel = m_waterGrid.at(Vector3i(i, j, k)).oldVelocity;
    double divergence = - vel[0] - vel[1] - vel[2];
    if (isInBoundsbyIdx(i + 1, j, k)) {
        divergence += m_waterGrid.at(Vector3i(i + 1, j, k)).oldVelocity[0];
    }
    if (isInBoundsbyIdx(i, j + 1, k)) {
        divergence += m_waterGrid.at(Vector3i(i, j + 1, k)).oldVelocity[1];
    }
    if (isInBoundsbyIdx(i, j, k + 1)) {
        divergence += m_waterGrid.at(Vector3i(i, j, k + 1)).oldVelocity[2];
    }
    return divergence;
}

Vector3d System::getCurl(int i, int j, int k){
    Vector3d curl(0, 0, 0);
    /// uz / y
    curl[0] += (j+1 < WATERGRID_Y) ? m_waterGrid[Vector3i(i, j+1, k)].oldVelocity[2] : 0;
    curl[0] -= (j-1 >= 0         ) ? m_waterGrid[Vector3i(i, j-1, k)].oldVelocity[2] : 0;
    /// uy / z
    curl[0] -= (k+1 < WATERGRID_X) ? m_waterGrid[Vector3i(i, j, k+1)].oldVelocity[1] : 0;
    curl[0] += (k-1 >= 0         ) ? m_waterGrid[Vector3i(i, j, k-1)].oldVelocity[1] : 0;

    /// ux / z
    curl[1] += (k+1 < WATERGRID_X) ? m_waterGrid[Vector3i(i, j, k+1)].oldVelocity[0] : 0;
    curl[1] -= (k-1 >= 0         ) ? m_waterGrid[Vector3i(i, j, k-1)].oldVelocity[0] : 0;
    /// uz / x
    curl[1] -= (i+1 < WATERGRID_X) ? m_waterGrid[Vector3i(i+1, j, k)].oldVelocity[2] : 0;
    curl[1] += (i-1 >= 0         ) ? m_waterGrid[Vector3i(i-1, j, k)].oldVelocity[2] : 0;

    /// uy / x
    curl[2] += (i+1 < WATERGRID_X) ? m_waterGrid[Vector3i(i+1, j, k)].oldVelocity[1] : 0;
    curl[2] -= (i-1 >= 0         ) ? m_waterGrid[Vector3i(i-1, j, k)].oldVelocity[1] : 0;
    /// ux / y
    curl[2] -= (j+1 < WATERGRID_Y) ? m_waterGrid[Vector3i(i, j+1, k)].oldVelocity[0] : 0;
    curl[2] += (j-1 >= 0         ) ? m_waterGrid[Vector3i(i, j-1, k)].oldVelocity[0] : 0;

    return curl;
}

Vector3d System::getCurlGradient(int i, int j, int k){
    Vector3d gradient(0, 0, 0);
    gradient[0] += (i+1 < WATERGRID_X) ? m_waterGrid[Vector3i(i+1, j, k)].curl.norm() : 0;
    gradient[0] -= (i-1 >= 0         ) ? m_waterGrid[Vector3i(i-1, j, k)].curl.norm() : 0;
    gradient[1] += (j+1 < WATERGRID_Y) ? m_waterGrid[Vector3i(i, j+1, k)].curl.norm() : 0;
    gradient[1] -= (j-1 >= 0         ) ? m_waterGrid[Vector3i(i, j-1, k)].curl.norm() : 0;
    gradient[2] += (k+1 < WATERGRID_X) ? m_waterGrid[Vector3i(i, j, k+1)].curl.norm() : 0;
    gradient[2] -= (k-1 >= 0         ) ? m_waterGrid[Vector3i(i, j, k-1)].curl.norm() : 0;
//    std::cout << gradient << std::endl;
    return gradient;
}

double isZero(double f) {
    return -EPSILON <= f && f <= EPSILON;
}
/**
 * Performs the Laplacian Operator on the velocity vector field.
 *
 * i, j, k : cell index
 * idx     : component of the velocity vector (either 0, 1, or 2 for x, y, or z)
 */
double System::laplacianOperatorOnVelocity(int i, int j, int k, int idx) {
    double laplacianVelocity = 0;

    /// i direction
    if (i + 1 < WATERGRID_X && !isZero(m_waterGrid.at(Vector3i{i+1, j, k}).oldVelocity[idx]))
        laplacianVelocity += m_waterGrid.at(Vector3i{i+1, j, k}).oldVelocity[idx];
    if (i - 1 >= 0 && !isZero(m_waterGrid.at(Vector3i{i-1, j, k}).oldVelocity[idx]))
        laplacianVelocity += m_waterGrid.at(Vector3i{i-1, j, k}).oldVelocity[idx];

    /// j direction
    if (j+1 < WATERGRID_Y &&  !isZero(m_waterGrid.at(Vector3i{i, j+1, k}).oldVelocity[idx]))
        laplacianVelocity += m_waterGrid.at(Vector3i{i, j+1, k}).oldVelocity[idx];
    if (j-1 >= 0 && !isZero(m_waterGrid.at(Vector3i{i, j-1, k}).oldVelocity[idx]))
        laplacianVelocity += m_waterGrid.at(Vector3i{i, j-1, k}).oldVelocity[idx];

    /// k direction
    if (k+1 < WATERGRID_Z && !isZero(m_waterGrid.at(Vector3i{i, j, k+1}).oldVelocity[idx]))
        laplacianVelocity += m_waterGrid.at(Vector3i{i, j, k+1}).oldVelocity[idx];
    if (k-1 >= 0 && !isZero(m_waterGrid.at(Vector3i{i, j, k-1}).oldVelocity[idx]))
        laplacianVelocity += m_waterGrid.at(Vector3i{i, j, k-1}).oldVelocity[idx];

    /// -6*currCellOldVelocity term
    if (!isZero(m_waterGrid.at(Vector3i{i, j, k}).oldVelocity[idx]))
        laplacianVelocity -= 6 * m_waterGrid.at(Vector3i{i, j, k}).oldVelocity[idx];

    return laplacianVelocity;
}


//// Get the interpolated velocity at a point in space.
Vector3d System::getVelocity(Vector3d pos){
    double x = getInterpolatedValue( pos.x() / CELL_DIM,         (pos.y() / CELL_DIM) - 0.5f, (pos.z() / CELL_DIM) - 0.5f, 0);
    double y = getInterpolatedValue((pos.x() / CELL_DIM) - 0.5f,  pos.y() / CELL_DIM,         (pos.z() / CELL_DIM) - 0.5f, 1);
    double z = getInterpolatedValue((pos.x() / CELL_DIM) - 0.5f, (pos.y() / CELL_DIM) - 0.5f,  pos.z() / CELL_DIM,         2);
    return Vector3d(x, y, z);
}

Eigen::Vector3d System::getVelocity(Eigen::Vector3d pos, CellBFECCField field) {
    double x = getInterpolatedValue( pos.x() / CELL_DIM,         (pos.y() / CELL_DIM) - 0.5f, (pos.z() / CELL_DIM) - 0.5f, 0, field);
    double y = getInterpolatedValue((pos.x() / CELL_DIM) - 0.5f,  pos.y() / CELL_DIM,         (pos.z() / CELL_DIM) - 0.5f, 1, field);
    double z = getInterpolatedValue((pos.x() / CELL_DIM) - 0.5f, (pos.y() / CELL_DIM) - 0.5f,  pos.z() / CELL_DIM,         2, field);
    return Vector3d(x, y, z);
}

std::vector<Vector3i> System::getGridNeighbors(int i, int j, int k){
    std::vector<Vector3i> neighbors;
    if (isInBoundsbyIdx(i + 1, j, k)) {
        neighbors.push_back(Vector3i(i + 1, j, k));
    } 
    if (isInBoundsbyIdx(i - 1, j, k)) {
        neighbors.push_back(Vector3i(i - 1, j, k));
    }
    if (isInBoundsbyIdx(i, j + 1, k)) {
        neighbors.push_back(Vector3i(i, j + 1, k));
    }
    if (isInBoundsbyIdx(i, j - 1, k)) {
        neighbors.push_back(Vector3i(i, j - 1, k));
    }
    if (isInBoundsbyIdx(i, j, k + 1)) {
        neighbors.push_back(Vector3i(i, j, k + 1));
    }
    if (isInBoundsbyIdx(i, j, k - 1)) {
        neighbors.push_back(Vector3i(i, j, k - 1));
    }
    return neighbors;
}

//// Checks whether a point position is within bounds of the waterGrid
bool System::isInBounds(double x, double y, double z) {
    // TODO: bounds checking here?
    bool xIs = (x >= 0) && (x < WATERGRID_X*CELL_DIM);
    bool yIs = (y >= 0) && (y < WATERGRID_Y*CELL_DIM);
    bool zIs = (z >= 0) && (z < WATERGRID_Z*CELL_DIM);
    return xIs && yIs && zIs;
}

//// Checks whether a cell index is within bounds of the waterGrid
bool System::isInBoundsbyIdx(int i, int j, int k) {
    bool iIs = (i >= 0) && (i < WATERGRID_X);
    bool jIs = (j >= 0) && (j < WATERGRID_Y);
    bool kIs = (k >= 0) && (k < WATERGRID_Z);
    return iIs && jIs && kIs;
}

Eigen::Vector3d System::getVelocityFromField(Eigen::Vector3i pos, CellBFECCField field) {
    switch (field) {
        case OLDVELOCITY:
            return this->m_waterGrid.at(pos).oldVelocity;
        case USTARFORWARD:
            return this->m_waterGrid.at(pos).uStarForward;
        case USTAR:
            return this->m_waterGrid.at(pos).uStar;
        case USQUIGGLY:
            return this->m_waterGrid.at(pos).uSquiggly;
        case CURRVELOCITY:
            return this->m_waterGrid.at(pos).currVelocity;
    }

    throw std::runtime_error("ERROR: CellBFECCField NOT FOUND");
}


//// Get an interpolated data value from the grid
double System::getInterpolatedValue(double x, double y, double z, int idx) {
    int i = floor(x);
    int j = floor(y);
    int k = floor(z);
    double weightAccum = 0;
    double totalAccum = 0;
    
    if (isInBoundsbyIdx(i, j, k)) {
        weightAccum += (i + 1 - x) * (j + 1 - y) * (k + 1 - z);
        totalAccum  += (i + 1 - x) * (j + 1 - y) * (k + 1 - z) * m_waterGrid.at(Vector3i(i, j, k)).oldVelocity[idx];
    }

    if (isInBoundsbyIdx(i + 1, j, k)) {
        totalAccum  += (x - i) * (j + 1 - y) * (k + 1 - z) * m_waterGrid.at(Vector3i(i + 1, j, k)).oldVelocity[idx];
        weightAccum += (x - i) * (j + 1 - y) * (k + 1 - z);
    }

    if (isInBoundsbyIdx(i, j+1, k)) {
        totalAccum  += (i + 1 - x) * (y - j) * (k + 1 - z) * m_waterGrid.at(Vector3i(i, j + 1, k)).oldVelocity[idx];
        weightAccum += (i + 1 - x) * (y - j) * (k + 1 - z);
    }

    if (isInBoundsbyIdx(i + 1, j+1, k)) {
        totalAccum  += (x - i) * (y - j) * (k + 1 - z) * m_waterGrid.at(Vector3i(i + 1, j + 1, k)).oldVelocity[idx];
        weightAccum += (x - i) * (y - j) * (k + 1 - z);
    }

    if (isInBoundsbyIdx(i, j, k+1)) {
        totalAccum  += (i + 1 - x) * (j + 1 - y) * (z - k) * m_waterGrid.at(Vector3i(i, j, k + 1)).oldVelocity[idx];
        weightAccum += (i + 1 - x) * (j + 1 - y) * (z - k);
    }

    if (isInBoundsbyIdx(i+1, j, k+1)) {
        totalAccum += (x - i) * (j + 1 - y) * (z - k) * m_waterGrid.at(Vector3i(i + 1, j, k + 1)).oldVelocity[idx];
        totalAccum += (x - i) * (j + 1 - y) * (z - k);
    }
    
    if (isInBoundsbyIdx(i, j+1, k+1)) {
        totalAccum  += (i + 1 - x) * (y - j) * (z - k) * m_waterGrid.at(Vector3i(i, j + 1, k + 1)).oldVelocity[idx];
        weightAccum += (i + 1 - x) * (y - j) * (z - k);
    }

    if (isInBoundsbyIdx(i+1, j+1, k+1)) {
        totalAccum  += (x - i) * (y - j) * (z - k) * m_waterGrid.at(Vector3i(i + 1, j + 1, k + 1)).oldVelocity[idx];
        weightAccum += (x - i) * (y - j) * (z - k);
    }

    if (weightAccum == 0) {
        return 0;
    }
    return totalAccum / weightAccum;
}

//// Get an interpolated data value from the grid with a givn field
double System::getInterpolatedValue(double x, double y, double z, int idx, CellBFECCField field) {
    int i = floor(x);
    int j = floor(y);
    int k = floor(z);
    double weightAccum = 0;
    double totalAccum = 0;

    if (isInBoundsbyIdx(i, j, k)) {
        weightAccum += (i + 1 - x) * (j + 1 - y) * (k + 1 - z);
        totalAccum  += (i + 1 - x) * (j + 1 - y) * (k + 1 - z) * getVelocityFromField(Vector3i(i, j, k), field)[idx];
    }

    if (isInBoundsbyIdx(i + 1, j, k)) {
        totalAccum  += (x - i) * (j + 1 - y) * (k + 1 - z) * getVelocityFromField(Vector3i(i + 1, j, k), field)[idx];
        weightAccum += (x - i) * (j + 1 - y) * (k + 1 - z);
    }

    if (isInBoundsbyIdx(i, j+1, k)) {
        totalAccum  += (i + 1 - x) * (y - j) * (k + 1 - z) * getVelocityFromField(Vector3i(i, j + 1, k), field)[idx];
        weightAccum += (i + 1 - x) * (y - j) * (k + 1 - z);
    }

    if (isInBoundsbyIdx(i + 1, j+1, k)) {
        totalAccum  += (x - i) * (y - j) * (k + 1 - z) * getVelocityFromField(Vector3i(i + 1, j + 1, k), field)[idx];
        weightAccum += (x - i) * (y - j) * (k + 1 - z);
    }

    if (isInBoundsbyIdx(i, j, k+1)) {
        totalAccum  += (i + 1 - x) * (j + 1 - y) * (z - k) * getVelocityFromField(Vector3i(i, j, k + 1), field)[idx];
        weightAccum += (i + 1 - x) * (j + 1 - y) * (z - k);
    }

    if (isInBoundsbyIdx(i+1, j, k+1)) {
        totalAccum += (x - i) * (j + 1 - y) * (z - k) * getVelocityFromField(Vector3i(i + 1, j, k + 1), field)[idx];
        totalAccum += (x - i) * (j + 1 - y) * (z - k);
    }

    if (isInBoundsbyIdx(i, j+1, k+1)) {
        totalAccum  += (i + 1 - x) * (y - j) * (z - k) * getVelocityFromField(Vector3i(i, j + 1, k + 1), field)[idx];
        weightAccum += (i + 1 - x) * (y - j) * (z - k);
    }

    if (isInBoundsbyIdx(i+1, j+1, k+1)) {
        totalAccum  += (x - i) * (y - j) * (z - k) * getVelocityFromField(Vector3i(i + 1, j + 1, k + 1), field)[idx];
        weightAccum += (x - i) * (y - j) * (z - k);
    }

    if (isZero(totalAccum) || isZero(weightAccum)) {
        return 0;
    }
    double total = totalAccum / weightAccum;
    return totalAccum / weightAccum;
}

/************************* DEBUGGING UTILS *****************************/
bool System::hasNan(Eigen::Vector3d v) {
    return isnan(v[0]) || isnan(v[1]) || isnan(v[2]);
}

bool System::hasInf(Eigen::Vector3d v) {
    return isinf(v[0]) || isinf(v[1]) || isinf(v[2]);
}

void System::checkNanAndInf() {
    // check watergrid
    #pragma omp parallel for
    for (auto& [k, v] : this->m_waterGrid) {
        auto k1 = k;
        auto v1 = v;
        assert(!hasNan(v.oldVelocity) && !hasInf(v.oldVelocity));
        assert(!hasNan(v.currVelocity) && !hasInf(v.currVelocity));
        assert(v.oldVelocity.norm() < 10000);
        assert(v.currVelocity.norm() < 10000);
    }

    // check ink
    #pragma omp parallel for
    for (auto& particle : this->m_ink) {
        assert(!hasNan(particle.position) && !hasInf(particle.position));
        assert(!hasNan(particle.velocity) && !hasInf(particle.position));
        assert(particle.velocity.norm() < 10000);
    }
}

/************************** PRINTING UTILS *****************************/

/// print a cell
ostream& operator<<(ostream& strm, const Cell& obj) {
    strm << "\tcurrent velocity: (" << obj.currVelocity.x() << ", ";
    strm << obj.currVelocity.y() << ", " << obj.currVelocity.z() << ")\n";
    strm << "\tcurl: <" << obj.curl.x() << "," << obj.curl.y() << "," << obj.curl.z() << ")\n";
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
        strm << "\tpos in hashmap: (" << k.x() << ", ";
        strm << k.x() << ", " << k.z() << ")\n";
        strm << v << endl;
    }

    strm << "********* PRINTING PARTICLES ***********\n";
    for (auto& el : obj.m_ink) {
        strm << el << endl;
    }
    return strm;
}
