#include <QCoreApplication>
#include <QCommandLineParser>

#include <iostream>
#include <QImage>
#include <fstream>


int main(int argc, char* argv[])
{
    // CLI Stuff
    QCoreApplication a(argc, argv);
    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addPositionalArgument("scene", "Scene file to be rendered");
    parser.addPositionalArgument("output", "Image file to write the rendered image to");
    parser.process(a);

    // TODO: Write to example file
    std::ofstream myfile;
    std::cout<<"writing this to a file"<<std::endl;
    myfile.open("output/EXAMPLEE.txt");
    myfile << "Writing this to a file. \n";
    myfile.close();


    a.exit();
}
