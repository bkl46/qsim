
#include "../include/qsim/qsim.hpp"
#include<pybind11/pybind11.h>



namespace py = pybind11;


PYBIND11_MODULE(qsim_py, m) {
    m.doc() = "Quantum Simulation Library - Path Integral Monte Carlo";

    // Register version
    m.attr("__version__") = "1.0.0";
    m.def("simulate", &qsim::simulate,
             "print out some simulation");

}







