#include "inkSim.h"
#include <fstream>
#include <chrono>

#define USE_BINARY false // binary or ascii

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
    // Get ink particles
    const std::vector<std::vector<Particle>> allInkParticles = this->ink_system.getInkParticles();
    #pragma omp parallel for
    for (int i = 0; i < allInkParticles.size(); i++) {
        // file for each drop
        std::string path = this->writeDirectory + "/drop" + std::to_string(i);
        std::filesystem::create_directory(path);
        std::string filename = path + "/" + std::to_string(frameNum) + ".ply";
        std::cout << "writing to: " + filename << std::endl;
        std::ofstream myfile;

        // write inks
        std::vector<Particle> inkParticles = allInkParticles[i];
        std::string numVertices = std::to_string(inkParticles.size());
        std::cout << "recording " + numVertices + " particles"<< std::endl;
        // setup header
        if (USE_BINARY) {
            myfile.open(filename, std::ios::binary | std::ios::out); // based off working directory not current file
            std::string header = "ply\n"
                                 "format binary_small_endian 1.0\n"
                                 "element vertex " + numVertices + "\n"
                                 "property float x\n"
                                 "property float y\n"
                                 "property float z\n"
                                 "end_header\n";
            myfile.write(header.c_str(), header.size());
        } else { /// ascii
            myfile.open(filename); // based off working directory not current file
            myfile << "ply\n";
            myfile << "format ascii 1.0\n";
            myfile << "element vertex " + numVertices + "\n";
            myfile << "property float x\n";
            myfile << "property float y\n";
            myfile << "property float z\n";
            myfile << "end_header\n";
        }

        // input data
        for(int i = 0; i<inkParticles.size(); i++) {
            Eigen::Vector3f currParticle = inkParticles[i].position;
            if (USE_BINARY) {
                myfile.write((char*)&currParticle.x(), sizeof(currParticle.x()));
                myfile.write((char*)&currParticle.y(), sizeof(currParticle.y()));
                myfile.write((char*)&currParticle.z(), sizeof(currParticle.z()));
            } else {
                Eigen::Vector3f currParticle = inkParticles[i].position;
                std::string position = std::to_string(currParticle.x()) + " " + std::to_string(currParticle.y()) + " " + std::to_string(currParticle.z());
                myfile << position + "\n";
            }
        }
        myfile.close();
    }
}

void InkSim::writeWaterGridVelocities() {
    std::string filename = this->writeDirectory + "/" + "velocityField.ply";
    std::ofstream myfile;

    std::cout << "writing to: " + filename << std::endl;
    std::string numVectors = std::to_string(WATERGRID_X * WATERGRID_Y * WATERGRID_Z);
    std::cout << "recording " + numVectors + " velocities" << std::endl;

    // set up header
    if (USE_BINARY) {
        myfile.open(filename, std::ios::binary | std::ios::out);
        std::string header = "ply\n"
                             "format binary_small_endian 1.0\n"
                             "element vertex " + numVectors + "\n"
                             "property float x\n"
                             "property float y\n"
                             "property float z\n"
                             "property float u\n"
                             "property float v\n"
                             "property float w\n"
                             "end_header\n";
        myfile.write(header.c_str(), header.size());
    } else { /// ascii
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
    }

    // Input velocity positions
    for (auto kv : ink_system.getWaterGrid()) {
        // Get cell center position
        float offset = CELL_DIM/2.f;
        Eigen::Vector3f cellCenter = Eigen::Vector3f{kv.first.x() + offset, kv.first.y() + offset, kv.first.z() + offset};

        // Get cell velocity
        Eigen::Vector3f cellVel = kv.second.currVelocity;

        if (USE_BINARY) {
            myfile.write((char*)&cellCenter.x(), sizeof(cellCenter.x()));
            myfile.write((char*)&cellCenter.y(), sizeof(cellCenter.y()));
            myfile.write((char*)&cellCenter.z(), sizeof(cellCenter.z()));

            myfile.write((char*)&cellVel.x(), sizeof(cellVel.x()));
            myfile.write((char*)&cellVel.y(), sizeof(cellVel.y()));
            myfile.write((char*)&cellVel.z(), sizeof(cellVel.z()));
        } else {
            // Convert cell position and velocity into a string
            std::string str = std::to_string(cellCenter.x()) + " " + std::to_string(cellCenter.y()) + " " + std::to_string(cellCenter.z()) + " "
                    + std::to_string(cellVel.x()) + " " + std::to_string(cellVel.y()) + " " + std::to_string(cellVel.z());

            // Write to file
            myfile << str + "\n";
        }
    }
    myfile.close();
}
