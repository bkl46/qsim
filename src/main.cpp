#include<iostream>
#include "../include/qsim/qsim.hpp"

using namespace qsim;


int main(int argc, char* argv[])
{

    std::cout<<"Hello World\n";

    if(argc > 1)
    {
        int i = atoi(argv[1]);
        std::cout << i << " entry \n";
    }



    simulate();

    XoshiroCpp::Xoshiro256StarStar rng(10);

    Particle p = Particle::boson(5);

    std::cout << "particle with \n";

    std::cout << "  " << p.mass() << " au mass\n";
    std::cout << "  " << p.charge() << " charge\n";
    std::cout << "  " << p.spin() << " spin\n";

    std::cout << "rng created, number: " << rng() << "\n";

    Bead<double> bead(1,2,3);

    std::cout << "bead created, at: " << bead[0] << ", " << bead[1] << ", " << bead[2] <<"\n";

    Path<double> path;

    std::cout << "path created, with " << path.num_beads() << " beads\n";

    potential::HarmonicOscillator harmonic(5);

    std::cout << "harmonic created\n";

    return 0;
}
