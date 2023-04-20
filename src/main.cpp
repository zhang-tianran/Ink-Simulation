#include <QCoreApplication>
#include <QCommandLineParser>

#include <iostream>
#include <QImage>
#include "inkSim.h"
#include <filesystem>


#define path "output" // Directory name to write to
int main(int argc, char* argv[])
{
    // CLI Stuff
    QCoreApplication a(argc, argv);
    QCommandLineParser parser;
    parser.addHelpOption();
    // parser.addPositionalArgument("numTimesteps", "Scene file to be rendered");
    // parser.addPositionalArgument("numTotalsteps", "Image file to write the rendered image to");
    parser.process(a);

    /// Create output directory if it doesn't already exist
    std::filesystem::create_directory(path);

    InkSim sim(path);
    sim.simulate(5, 10); // TODO: input stuff potentially
    a.exit();
}
