#include "inkSim.h"
#include <fstream>
#include <chrono>

InkSim::InkSim(std::string writeDirectory) {
    this->ink_system.init();
    this->ink_system.checkNanAndInf();
    this->writeDirectory = writeDirectory;
}

void InkSim::simulate(const float renderTimestep, int totalTimesteps) {
    /// For recording time
    auto startTS = std::chrono::system_clock::now();
    auto startTime = std::chrono::system_clock::to_time_t(startTS);
    std::cout << "\e[32mSimulation started at: \e[m" << std::ctime(&startTime) << std::endl;

    int timestepCounter = 1;
    int currFrame = 1;
    float timeSinceLastRender = 0;
    while (timestepCounter <= totalTimesteps) {
        float timeStep = this->ink_system.solve(renderTimestep - timeSinceLastRender);
//        std::cout << "TIMESTEP INK: " << timeStep << std::endl;
//        std::cout << "NUM TIMESTEPS: " << timestepCounter << std::endl;
        timeSinceLastRender += timeStep;
        if (timeSinceLastRender>=renderTimestep) {
            writeToFile(currFrame);
            currFrame+=1;
            timeSinceLastRender = 0;
        }
        timestepCounter += 1;
    }

    writeWaterGridVelocities();

    /// For recording time
    auto endTS = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endTS-startTS;
    int seconds = elapsed_seconds.count();
    std::cout << "\e[32mSimulation & Writing took " << seconds / 60 << " min " << seconds % 60 << " sec\e[m" << std::endl;
}


void InkSim::writeToFile(int frameNum) {
    std::string filename = this->writeDirectory + "/" + std::to_string(frameNum) + ".ply";
    std::ofstream myfile;
    std::cout << "writing to: " + filename << std::endl;
    // Get ink particles
    const std::vector<Particle> inkParticles = this->ink_system.getInkParticles();
    std::string numVertices = std::to_string(inkParticles.size());
    std::cout << "recording " + numVertices + " particles"<< std::endl;
    myfile.open(filename); // based off working directory not current file
    // setup header
    myfile << "ply\n";
    myfile << "format ascii 1.0\n";
    myfile << "element vertex " + numVertices + "\n";
    myfile << "property float x\n";
    myfile << "property float y\n";
    myfile << "property float z\n";
    myfile << "end_header\n";
    // input data
    for(int i = 0; i<inkParticles.size(); i++) {
        Eigen::Vector3f currParticle = inkParticles[i].position;
        std::string position = std::to_string(currParticle.x()) + " " + std::to_string(currParticle.y()) + " " + std::to_string(currParticle.z());
        myfile << position + "\n";
    }
    myfile.close();
}

void InkSim::writeWaterGridVelocities() {
    std::string filename = this->writeDirectory + "/" + "velocityField.ply";
    std::ofstream myfile;

    std::cout << "writing to: " + filename << std::endl;
    std::string numVectors = std::to_string(WATERGRID_X * WATERGRID_Y * WATERGRID_Z);
    std::cout << "recording " + numVectors + " velocities" << std::endl;

    // Set up header for velocityPositions.ply
    myfile.open(filename);
    myfile << "ply\n";
    myfile << "format ascii 1.0\n";
    myfile << "element vertex " + numVectors + "\n";
    myfile << "property float x\n";
    myfile << "property float y\n";
    myfile << "property float z\n";
    myfile << "property float u\n";
    myfile << "property float v\n";
    myfile << "property float w\n";
    myfile << "end_header\n";

    // Input velocity positions
    for (auto kv : ink_system.getWaterGrid()) {
        // Get cell center position
        float offset = CELL_DIM/2.f;
        Eigen::Vector3f cellCenter = Eigen::Vector3f{kv.first.x() + offset, kv.first.y() + offset, kv.first.z() + offset};

        // Get cell velocity
        Eigen::Vector3f cellVel = kv.second.currVelocity;

        // Convert cell position and velocity into a string
        std::string str = std::to_string(cellCenter.x()) + " " + std::to_string(cellCenter.y()) + " " + std::to_string(cellCenter.z()) + " "
                + std::to_string(cellVel.x()) + " " + std::to_string(cellVel.y()) + " " + std::to_string(cellVel.z());

        // Write to file
        myfile << str + "\n";
    }
    myfile.close();
}


void InkSim::writeToFile() {
    // TODO: Write to actual file with actualy data
    std::string filename = this->writeDirectory + "/example.ply";
    std::ofstream myfile;
    std::cout << "writing to: " + filename << std::endl;
    myfile.open(filename); // based off working directory not current file
    // setup header
    myfile << "ply\n";
    myfile << "format ascii 1.0\n";
    myfile << "element vertex 3\n";
    myfile << "property float x\n";
    myfile << "property float y\n";
    myfile << "property float z\n";
    myfile << "end_header\n";
    // input data
    myfile << "1 0 0\n";
    myfile << "0 1 0\n";
    myfile << "0 0 1\n";
    myfile.close();
}
