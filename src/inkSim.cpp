#include "inkSim.h"

InkSim::InkSim(std::string fileDirectory) {
    this->ink_system = System();
    this->ink_system.init();
    this->fileDirectory = fileDirectory;
}

void simulate(int numTimesteps, int totalTimesteps) {
    // TODO
}

void writeToFile() {
    // TODO: Write to actual file with actualy data
    std::string filename = this->fileDirectory + "/example.ply";
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