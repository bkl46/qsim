// src/core/types_demo.cpp (optional, for documentation)
#include "../include/qsim/core/types.hpp"
#include <iostream>

// Example of using the strong types
void demonstrate_types() {
    using namespace qsim;
    
    // Strong types prevent accidental mixing
    Mass m1(1.0);
    Mass m2(2.0);
    Length l1(3.0);
    
    // This would compile but is semantically wrong:
    // auto wrong = m1 + l1;  // ERROR: different types cannot be added
    
    // Correct operations preserve types
    auto m3 = m1 + m2;  // Mass
    auto m4 = m1 * 2.0; // Mass
    
    // Vector operations
    Vec3d pos(1.0, 2.0, 3.0);
    Vec3d vel(0.1, 0.2, 0.3);
    auto momentum = pos * 0.5;  // Scaling
    
    
    // Simulation parameters
    SimulationParameters params(1.0, 100, 1000);
    std::cout << "Time step: " << params.time_step << std::endl;
}


int main(){
    demonstrate_types();
    return 0;
}
