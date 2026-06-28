# tests/CMakeLists.txt

enable_testing()

set(TEST_SOURCES
    cpp/unit/test_types.cpp
    cpp/unit/test_parallel.cpp
)

add_executable(qsim_tests ${TEST_SOURCES})

target_link_libraries(qsim_tests
    PUBLIC
    Catch2::Catch2WithMain
    qsim_core 
)

#target_include_directories(qsim_tests
#PRIVATE
#${CMAKE_SOURCE_DIR}/include
#${CMAKE_SOURCE_DIR}/src
#)

# register test
add_test(NAME qsim_tests COMMAND qsim_tests)




add_custom_target(check
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    DEPENDS qsim_tests
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running all tests with CTest"
)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16)
    include(Catch)
    catch_discover_tests(qsim_tests)
endif()









