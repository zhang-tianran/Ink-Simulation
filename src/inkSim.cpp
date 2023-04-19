#include "inkSim.h"
#include <fstream>

InkSim::InkSim(std::string writeDirectory) {
    this->ink_system = System();
    this->ink_system.init();
    this->writeDirectory = writeDirectory;
}

void InkSim::simulate(int numTimesteps, int totalTimesteps) {
    // TODO
    int timestepCounter = 1;
    while (timestepCounter <= totalTimesteps) {
        if (timestepCounter % numTimesteps == 0) {
            writeToFile(timestepCounter/numTimesteps);
        }
        timestepCounter += 1;
    }

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
