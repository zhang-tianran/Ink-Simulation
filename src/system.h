#pragma once
#include <iostream>
#include <unordered_map>
#include <Eigen/Dense>

typedef struct Cell {
    Eigen::Vector3f velocity;
    float pressure;
} Cell;

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

    std::unordered_map<Eigen::Vector3f, Cell, hash_func> waterGrid;
    float updateWaterGrid(std::unordered_map<Eigen::Vector3f, Cell, hash_func> &waterGrid);
    float calcTimeStep(const std::unordered_map<Eigen::Vector3f, Cell, hash_func> &waterGrid);
    void updateGridFromMarkers(std::unordered_map<Eigen::Vector3f, Cell, hash_func> &waterGrid);
    void updateVelocityField(std::unordered_map<Eigen::Vector3f, Cell, hash_func> &waterGrid);

};


