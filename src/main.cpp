#include<iostream>
#include "../include/qsim/qsim.hpp"

using namespace qsim;


int main(int argc, char* argv[])
{

    std::cout<<"Hello World\n";



    simulate();

    XoshiroCpp::Xoshiro256StarStar rng(10);

    std::cout << "rng created, number: " << rng() << "\n";

    Bead<double> bead(1,2,3);

    std::cout << "bead created, at: " << bead[0] << ", " << bead[1] << ", " << bead[2] <<"\n";

    Path<double> path;

    std::cout << "path created, with " << path.num_beads() << " beads\n";

    potential::HarmonicOscillator harmonic(5);

    std::cout << "harmonic created\n";

    return 0;
}
