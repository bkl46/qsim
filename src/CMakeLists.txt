
set(QSIM_SOURCES
    core/simulator.cpp
    main.cpp
)

add_library(qsim_core  STATIC
    ${QSIM_SOURCES}
)


set_target_properties(qsim_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(qsim_core PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)


add_executable(qsim_cli main.cpp)

target_link_libraries(qsim_cli PRIVATE qsim_core)



