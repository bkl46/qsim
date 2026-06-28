
if(pybind11_FOUND)

    pybind11_add_module(qsim_py
        qsim_py.cpp
    )

    target_link_libraries(qsim_py PUBLIC qsim_core)


    install(
        TARGETS qsim_py
        LIBRARY DESTINATION qsim   
        RUNTIME DESTINATION qsim    
    )


endif()

