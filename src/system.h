#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>

typedef struct Cell {
    Eigen::Vector3f velocity;
    float pressure;
} Cell;

typedef struct Particle {
    float opacity;
    float lifeTime;
    Eigen::Vector3f velocity;
    Eigen::Vector3f position;
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

    std::unordered_map<Eigen::Vector3f, Cell, hash_func> waterGrid;
    float updateWaterGrid();
    float calcTimeStep();
    void updateGridFromMarkers();
    void updateVelocityField();

    std::vector<Particle> ink;
    void updateParticles(float timeStep);

};


