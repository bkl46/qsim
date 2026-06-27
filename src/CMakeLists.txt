#collect source files (you can list them explicitly or use GLOB)


set(QSIM_SOURCES
    core/simulator.cpp
    main.cpp
)

add_library(qsim_core  SHARED
    ${QSIM_SOURCES}
)

target_include_directories(qsim_core PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)


add_executable(qsim_cli main.cpp)

target_link_libraries(qsim_cli PRIVATE qsim_core)



