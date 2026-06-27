# Optional: built only if pybind11 is found
if(pybind11_FOUND)

    pybind11_add_module(qsim_py
        qsim_py.cpp
    )

    target_link_libraries(qsim_py PUBLIC qsim_core)

endif()
