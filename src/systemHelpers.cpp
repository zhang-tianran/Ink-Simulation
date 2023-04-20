#include "system.h"
using namespace Eigen;
using namespace std;

Vector3f System::getGradient(int i, int j, int k, VectorXf g){
    int row_idx = grid2mat(i, j, k);
    float val = g[row_idx];
    Vector3f gradient(val, val, val);
    if (isInBoundsbyIdx(i - 1, j, k)) {
        gradient[0] -= g[grid2mat(i - 1, j, k)];
    }
    if (isInBoundsbyIdx(i, j - 1, k)) {
        gradient[0] -= g[grid2mat(i, j - 1, k)];
    }
    if (isInBoundsbyIdx(i, j, k - 1)) {
        gradient[2] -= g[grid2mat(i, j, k - 1)];
    }
    return gradient;
}

float System::getDivergence(int i, int j, int k){
    Vector3f vel = m_waterGrid[Vector3i(i, j, k)].currVelocity;
    float divergence = - vel[0] - vel[1] - vel[2];
    if (isInBoundsbyIdx(i + 1, j, k)) {
        divergence += m_waterGrid[Vector3i(i + 1, j, k)].currVelocity[0];
    }
    if (isInBoundsbyIdx(i, j + 1, k)) {
        divergence += m_waterGrid[Vector3i(i, j + 1, k)].currVelocity[1];
    }
    if (isInBoundsbyIdx(i, j, k + 1)) {
        divergence += m_waterGrid[Vector3i(i, j, k + 1)].currVelocity[2];
    }
    return divergence;
}

//// Get the interpolated velocity at a point in space.
Vector3f System::getVelocity(Vector3f pos){
    float x = getInterpolatedValue(pos[0] / CELL_DIM,        pos[1] / CELL_DIM - 0.5f, pos[2] / CELL_DIM - 0.5f, 0);
    float y = getInterpolatedValue(pos[0] / CELL_DIM - 0.5f, pos[1] / CELL_DIM,        pos[2] / CELL_DIM - 0.5f, 1);
    float z = getInterpolatedValue(pos[0] / CELL_DIM - 0.5f, pos[1] / CELL_DIM - 0.5f, pos[2] / CELL_DIM,        2);
    return Vector3f(x, y, z);
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
bool System::isInBounds(float x, float y, float z) {
    // TODO: bounds checking here?
    bool xIs = (x >= 0) && (x <= WATERGRID_X*CELL_DIM);
    bool yIs = (y >= 0) && (y <= WATERGRID_Y*CELL_DIM);
    bool zIs = (z >= 0) && (z <= WATERGRID_Z*CELL_DIM);
    return xIs && yIs && zIs;
}

//// Checks whether a cell index is within bounds of the waterGrid
bool System::isInBoundsbyIdx(int i, int j, int k) {
    bool iIs = (i >= 0) && (i < WATERGRID_X);
    bool jIs = (j >= 0) && (j < WATERGRID_Y);
    bool kIs = (k >= 0) && (k < WATERGRID_Z);
    return iIs && jIs && kIs;
}

//// Get an interpolated data value from the grid
float System::getInterpolatedValue(float x, float y, float z, int idx) {
    int i = floor(x);
    int j = floor(y);
    int k = floor(z);
    float weightAccum = 0;
    float totalAccum = 0;
    
    if (isInBounds(i, j, k)) {
        weightAccum += (i + 1 - x) * (j + 1 - y) * (k + 1 - z);
        totalAccum  += (i + 1 - x) * (j + 1 - y) * (k + 1 - z) * m_waterGrid[Vector3i(i, j, k)].oldVelocity[idx];
    }

    if (isInBounds(i + 1, j, k)) {
        totalAccum  += (x - i) * (j + 1 - y) * (k + 1 - z) * m_waterGrid[Vector3i(i + 1, j, k)].oldVelocity[idx];
        weightAccum += (x - i) * (j + 1 - y) * (k + 1 - z);
    }

    if (isInBounds(i, j+1, k)) {
        totalAccum  += (i + 1 - x) * (y - j) * (k + 1 - z) * m_waterGrid[Vector3i(i, j + 1, k)].oldVelocity[idx];
        weightAccum += (i + 1 - x) * (y - j) * (k + 1 - z);
    }

    if (isInBounds(i + 1, j+1, k)) {
        totalAccum  += (x - i) * (y - j) * (k + 1 - z) * m_waterGrid[Vector3i(i + 1, j + 1, k)].oldVelocity[idx];
        weightAccum += (x - i) * (y - j) * (k + 1 - z);
    }

    if (isInBounds(i, j, k+1)) {
        totalAccum  += (i + 1 - x) * (j + 1 - y) * (z - k) * m_waterGrid[Vector3i(i, j, k + 1)].oldVelocity[idx];
        weightAccum += (i + 1 - x) * (j + 1 - y) * (z - k);
    }

    if (isInBounds(i+1, j, k+1)) {
        totalAccum += (x - i) * (j + 1 - y) * (z - k) * m_waterGrid[Vector3i(i + 1, j, k + 1)].oldVelocity[idx];
        totalAccum += (x - i) * (j + 1 - y) * (z - k);
    }
    
    if (isInBounds(i, j+1, k+1)) {
        totalAccum  += (i + 1 - x) * (y - j) * (z - k) * m_waterGrid[Vector3i(i, j + 1, k + 1)].oldVelocity[idx];
        weightAccum += (i + 1 - x) * (y - j) * (z - k);
    }

    if (isInBounds(i+1, j+1, k+1)) {
        totalAccum  += (x - i) * (y - j) * (z - k) * m_waterGrid[Vector3i(i + 1, j + 1, k + 1)].oldVelocity[idx];
        weightAccum += (x - i) * (y - j) * (z - k);
    }

    return totalAccum / weightAccum;
}
