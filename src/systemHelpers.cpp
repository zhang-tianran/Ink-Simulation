#include "system.h"
using namespace Eigen;
using namespace std;

Vector3f System::getGradient(int l, int w, int h, MatrixXf g){
    int row_idx = grid2mat(l, w, h);
    int val = g.row(row_idx);
    Vector3f gradient(val, val, val);
    if (isInBoundsbyIdx(l - 1, w, h)) {
        gradient[0] -= g.row(grid2mat(l - 1, w, h));
    }
    if (isInBoundsbyIdx(l, w - 1, h)) {
        gradient[0] -= g.row(grid2mat(l, w - 1, h));
    }
    if (isInBoundsbyIdx(l, w, h - 1)) {
        gradient[2] -= g.row(grid2mat(l, w, h - 1));
    }
    return gradient;
}

float System::getDivergence(int l, int w, int h){
    Vector3f vel = m_waterGrid[Vector3i(l, w, h)].velocity;
    float divergence = - vel[0] - vel[1] - vel[2];
    if (isInBoundsbyIdx(l + 1, w, h)) {
        divergence += m_waterGrid[Vector3i(l + 1, w, h)].velocity[0];
    }
    if (isInBoundsbyIdx(l, w + 1, h)) {
        divergence += m_waterGrid[Vector3i(l, w + 1, h)].velocity[1];
    }
    if (isInBoundsbyIdx(l, w, h + 1)) {
        divergence += m_waterGrid[Vector3i(l, w, h + 1)].velocity[2];
    }
    return divergence;
}

//// Get the interpolated velocity at a point in space.
Vector3f System::getVelocity(Vector3f pos){
    float x = getInterpolatedValue(pos[0] / CELL_DIM, pos[1] / CELL_DIM - 0.5f, pos[2] / CELL_DIM - 0.5f, 0);
    float y = getInterpolatedValue(pos[0] / CELL_DIM - 0.5f, pos[1] / CELL_DIM, pos[2] / CELL_DIM - 0.5f, 1);
    float z = getInterpolatedValue(pos[0] / CELL_DIM - 0.5f, pos[1] / CELL_DIM - 0.5f, pos[2] / CELL_DIM, 2);
    return Vector3f(x, y, z);
}

std::vector<Vector3i> System::getGridNeighbors(int l, int w, int h){
    std::vector<Vector3i> neighbors;
    if (isInBoundsbyIdx(l + 1, w, h)) {
        neighbors.push_back(Vector3i(l + 1, w, h));
    } 
    if (isInBoundsbyIdx(l - 1, w, h)) {
        neighbors.push_back(Vector3i(l - 1, w, h));                    
    }
    if (isInBoundsbyIdx(l, w + 1, h)) {
        neighbors.push_back(Vector3i(l, w + 1, h));
    }
    if (isInBoundsbyIdx(l, w - 1, h)) {
        neighbors.push_back(Vector3i(l, w - 1, h));   
    }
    if (isInBoundsbyIdx(l, w, h + 1)) {
        neighbors.push_back(Vector3i(l, w, h + 1));
    }
    if (isInBoundsbyIdx(l, w, h - 1)) {
        neighbors.push_back(Vector3i(l, w, h - 1));
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
bool System::isInBoundsbyIdx(int l, int w, int h) {
    bool xIs = (l >= 0) && (l < WATERGRID_X);
    bool yIs = (w >= 0) && (w < WATERGRID_Y);
    bool zIs = (h >= 0) && (h < WATERGRID_Z);
    return xIs && yIs && zIs;
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
