#ifndef INK_SIM_H
#define INK_SIM_H
#include "system.h"

class InkSim {
public:

    /**
     * @param fileDirectory, the path to the *directory* to save ply files to
    */
    InkSim(std::string writeDirectory);

    /**
     * run the ink simulation
     *
     * @param numTimesteps, the number of timesteps to do before writing to a file
     * @param totalTimesteps, total number of timesteps to simulate
    */
    void simulate(int numTimesteps, int totalTimesteps);

private:
    System ink_system;
    std::string writeDirectory;

    void writeToFile();
    void writeToFile(int frameNum);

    void writeWaterGridVelocities();
};

#endif // INK_SIM_H
