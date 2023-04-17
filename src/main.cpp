#include <QCoreApplication>
#include <QCommandLineParser>

#include <iostream>
#include <QImage>
#include "inkSim.h"


#define path "TODO: input something meaningful here"
int main(int argc, char* argv[])
{
    // CLI Stuff
    QCoreApplication a(argc, argv);
    QCommandLineParser parser;
    parser.addHelpOption();
    // parser.addPositionalArgument("numTimesteps", "Scene file to be rendered");
    // parser.addPositionalArgument("numTotalsteps", "Image file to write the rendered image to");
    parser.process(a);

    InkSim sim(path);
    sim.simulate(0, 0); // TODO: input stuff potentially
    a.exit();
}
