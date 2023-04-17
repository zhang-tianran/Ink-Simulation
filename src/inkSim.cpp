#include "inkSim.h"
#include <fstream>

InkSim::InkSim(std::string writeDirectory) {
    this->ink_system = System();
    this->ink_system.init();
    this->writeDirectory = writeDirectory;
}

void InkSim::simulate(int numTimesteps, int totalTimesteps) {
    // TODO
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
