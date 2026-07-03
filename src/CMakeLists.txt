
set(QSIM_SOURCES
    core/simulator.cpp
    main.cpp
)

add_library(qsim_core  STATIC
    ${QSIM_SOURCES}
)


set_target_properties(qsim_core PROPERTIES POSITION_INDEPENDENT_CODE ON)


if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(qsim_core PRIVATE
        $<$<CONFIG:Debug>:-O0 -g>
        $<$<CONFIG:Release>:-O3 -march=native -mtune=native>
        $<$<CONFIG:RelWithDebInfo>:-O3 -march=native -mtune=native -g>
    )
endif()

target_include_directories(qsim_core PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)


add_executable(qsim_cli main.cpp)

target_link_libraries(qsim_cli PRIVATE qsim_core)



