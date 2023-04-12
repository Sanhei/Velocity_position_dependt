#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include "options.hpp"
#include "get_rid_of_transition.hpp"

#define PROGRAM_NAME "D_trans"
int main(int argc, char** argv)
{
        // Read trajectory
        program_options options;
        options.parse(argc, argv);

        std::cout<<"File name is "<< options.get("filename")<<std::endl;
        std::string filename=  options.get("filename");//= options.get("filename");
        std::ifstream inputfile;

        inputfile.open(filename);
        std::cout<<"Check point One"<<std::endl;
        if(!inputfile)
        {
                std::cout<<"unable to openfile";
                exit(1);
        }
        double x;
        std::vector<double> trajectory;
        std::cout<<"Reading trajectory"<<std::endl;
        while(inputfile>>x)
        {
                trajectory.push_back(x);
        }

        double timestep = 0.001;
        int Relax_time = 10;
        int N_position = 100;
        int N_velocity = 100;
        GetRidofTransition D_trans(timestep, Relax_time, trajectory, N_position, N_velocity);
        std::cout<<"Transition calculation begins"<<std::endl;
        D_trans.Set_bins();
        std::cout<<"Deleting transition"<<std::endl;
        D_trans.Locate_transition();
        std::string save_file="./heatmap.txt";
        std::cout<<"Saving file"<<std::endl;
        D_trans.Save_data(save_file);
        return 0;
}


void program_options::parse(int argc, char** argv)
{
    namespace po = boost::program_options;

    po::options_description desc("Program options");
    desc.add_options()
        ("filename", po::value<std::string>()->required(), "trajectory filename")

        ("help", "display this help and exit")
        ;

    try {
        po::command_line_parser parser(argc, argv);
        po::parsed_options parsed(parser.options(desc).run());
        po::store(parsed, vm_);
        

        if (vm_.count("help")) {
            std::cout << "Usage: " PROGRAM_NAME " [OPTION]...\n\n" << desc << std::endl;
            exit(EXIT_SUCCESS);
        }

        po::notify(vm_);
    }
    catch (std::exception const& e) {
        std::cerr << PROGRAM_NAME ": " << e.what() << "\n\n";
        std::cerr << "Try `" PROGRAM_NAME " --help' for more information." << std::endl;
        exit(EXIT_FAILURE);
    }
}

