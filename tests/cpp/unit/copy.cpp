// test/cpp/types_test.cpp
#include "../../../include/qsim/core/types.hpp"
#include <iostream>
#include <cassert>
#include <type_traits>

using namespace qsim;

//helper
void print_test_result(const std::string& test_name, bool passed) {
    std::cout << "[ " << (passed ? "✓" : "✗") << " ] " << test_name << std::endl;
}

// Test 1: Strong Type Basic Operations
void test_strong_type_basic() {
    std::cout << "\n=== Testing Strong Type Basic Operations ===" << std::endl;
    
    // Construction
    Mass m1{1.0};
    Mass m2{2.0};
    Mass m3{3.0};
    
    // Addition
    Mass m_sum = m1 + m2;
    assert(m_sum.value() == 3.0);
    print_test_result("Strong type addition", m_sum.value() == 3.0);
    
    // Subtraction
    Mass m_diff = m3 - m1;
    assert(m_diff.value() == 2.0);
    print_test_result("Strong type subtraction", m_diff.value() == 2.0);
    
    // Multiplication
    Mass m_mul = m1 * m2;
    assert(m_mul.value() == 2.0);
    print_test_result("Strong type multiplication", m_mul.value() == 2.0);
    
    // Division
    Mass m_div = m3 / m1;
    assert(m_div.value() == 3.0);
    print_test_result("Strong type division", m_div.value() == 3.0);
    
    // Scalar multiplication (left)
    Mass m_scalar_left = 2.0 * m1;
    assert(m_scalar_left.value() == 2.0);
    print_test_result("Scalar multiplication (left)", m_scalar_left.value() == 2.0);
    
    // Scalar multiplication (right)
    Mass m_scalar_right = m1 * 2.0;
    assert(m_scalar_right.value() == 2.0);
    print_test_result("Scalar multiplication (right)", m_scalar_right.value() == 2.0);
    
    // Implicit conversion to underlying type
    double raw_value = m1;  // Should work due to operator T()
    assert(raw_value == 1.0);
    print_test_result("Implicit conversion to underlying type", raw_value == 1.0);
    
    // Value accessor
    assert(m1.value() == 1.0);
    print_test_result("Value accessor", m1.value() == 1.0);
    
    // Default construction
    Mass m_default;
    assert(m_default.value() == 0.0);
    print_test_result("Default construction", m_default.value() == 0.0);
}

// Test 2: Type Safety
void test_type_safety() {
    std::cout << "\n=== Testing Type Safety ===" << std::endl;
    
    Mass mass{5.0};
    Length length{3.0};
    Time time{2.0};
    
    // These should compile and work (same type operations)
    Mass mass_result = mass + mass;
    assert(mass_result.value() == 10.0);
    print_test_result("Same type addition compiles", true);
    
    // This should NOT compile if you uncomment it (type safety)
    // auto invalid = mass + length;  // This should cause compilation error
    
    // This should compile - scalar multiplication works
    auto scaled_mass = mass * 2.0;
    assert(scaled_mass.value() == 10.0);
    print_test_result("Scalar multiplication with different types", true);
    
    // Test that different strong types are indeed different types
    bool same_type = std::is_same_v<Mass, Length>;
    print_test_result("Different strong types are distinct types", !same_type);
    
    // Test StrongTypeConcept
    bool satisfies_concept = StrongTypeConcept<Mass>;
    print_test_result("StrongTypeConcept satisfied", satisfies_concept);
}

// Test 3: Vec3 Operations
void test_vec3_operations() {
    std::cout << "\n=== Testing Vec3 Operations ===" << std::endl;
    
    // Construction
    Vec3d v1{1.0, 2.0, 3.0};
    Vec3d v2{4.0, 5.0, 6.0};
    Vec3d v_scalar{2.0};  // Explicit constructor
    
    // Addition
    Vec3d v_sum = v1 + v2;
    assert(v_sum.x == 5.0 && v_sum.y == 7.0 && v_sum.z == 9.0);
    print_test_result("Vec3 addition", true);
    
    // Subtraction
    Vec3d v_diff = v2 - v1;
    assert(v_diff.x == 3.0 && v_diff.y == 3.0 && v_diff.z == 3.0);
    print_test_result("Vec3 subtraction", true);
    
    // Scalar multiplication
    Vec3d v_scaled = v1 * 2.0;
    assert(v_scaled.x == 2.0 && v_scaled.y == 4.0 && v_scaled.z == 6.0);
    print_test_result("Vec3 scalar multiplication", true);
    
    // Scalar division
    Vec3d v_div = v2 / 2.0;
    assert(v_div.x == 2.0 && v_div.y == 2.5 && v_div.z == 3.0);
    print_test_result("Vec3 scalar division", true);
    
    // In-place operators
    Vec3d v_plus_eq = v1;
    v_plus_eq += v2;
    assert(v_plus_eq.x == 5.0 && v_plus_eq.y == 7.0 && v_plus_eq.z == 9.0);
    print_test_result("Vec3 += operator", true);
    
    Vec3d v_minus_eq = v2;
    v_minus_eq -= v1;
    assert(v_minus_eq.x == 3.0 && v_minus_eq.y == 3.0 && v_minus_eq.z == 3.0);
    print_test_result("Vec3 -= operator", true);
    
    // Dot product
    double dot_result = v1.dot(v2);
    assert(dot_result == 1.0*4.0 + 2.0*5.0 + 3.0*6.0);  // 4 + 10 + 18 = 32
    print_test_result("Vec3 dot product", dot_result == 32.0);
    
    // Cross product
    Vec3d v_cross = v1.cross(v2);
    // v1 x v2 = (2*6 - 3*5, 3*4 - 1*6, 1*5 - 2*4) = (12-15, 12-6, 5-8) = (-3, 6, -3)
    assert(v_cross.x == -3.0 && v_cross.y == 6.0 && v_cross.z == -3.0);
    print_test_result("Vec3 cross product", true);
    
    // Norm squared
    Vec3d v{3.0, 4.0, 0.0};
    double norm_sq = v.norm_sq();
    assert(norm_sq == 25.0);
    print_test_result("Vec3 norm squared", norm_sq == 25.0);
    
    // Norm
    double norm = v.norm();
    assert(norm == 5.0);
    print_test_result("Vec3 norm", norm == 5.0);
    
    // Accessor
    assert(v[0] == 3.0 && v[1] == 4.0 && v[2] == 0.0);
    print_test_result("Vec3 operator[] access", true);
    
    // Equality
    Vec3d v_equal{3.0, 4.0, 0.0};
    assert(v == v_equal);
    print_test_result("Vec3 equality operator", v == v_equal);
    
    // Vector type aliases
    Vec3f v_float{1.0f, 2.0f, 3.0f};
    bool alias_works = std::is_same_v<decltype(v_float), Vec3<float>>;
    print_test_result("Vec3 type aliases (Vec3f)", alias_works);
}

// Test 4: ParticleIndex and BeadIndex
void test_indices() {
    std::cout << "\n=== Testing ParticleIndex and BeadIndex ===" << std::endl;
    
    // ParticleIndex
    ParticleIndex p1{5};
    ParticleIndex p2{5};
    ParticleIndex p3{10};
    
    assert(static_cast<i32>(p1) == 5);
    print_test_result("ParticleIndex to i32 conversion", static_cast<i32>(p1) == 5);
    
    assert(static_cast<usize>(p1) == 5);
    print_test_result("ParticleIndex to usize conversion", static_cast<usize>(p1) == 5);
    
    assert(p1 == p2);
    assert(p1 != p3);
    print_test_result("ParticleIndex comparison", p1 == p2 && p1 != p3);
    
    assert(p1 != nullptr);
    print_test_result("ParticleIndex nullptr comparison (valid)", p1 != nullptr);
    
    ParticleIndex p_null;
    assert(p_null == nullptr);
    print_test_result("ParticleIndex nullptr comparison (null)", p_null == nullptr);
    
    // Increment
    ParticleIndex p_inc{0};
    ++p_inc;
    assert(static_cast<i32>(p_inc) == 1);
    print_test_result("ParticleIndex prefix increment", static_cast<i32>(p_inc) == 1);
    
    ParticleIndex p_post{0};
    p_post++;
    assert(static_cast<i32>(p_post) == 1);
    print_test_result("ParticleIndex postfix increment", static_cast<i32>(p_post) == 1);
    
    // BeadIndex
    BeadIndex b1{5};
    BeadIndex b2{5};
    BeadIndex b3{10};
    
    assert(static_cast<i32>(b1) == 5);
    print_test_result("BeadIndex to i32 conversion", static_cast<i32>(b1) == 5);
    
    assert(static_cast<usize>(b1) == 5);
    print_test_result("BeadIndex to usize conversion", static_cast<usize>(b1) == 5);
    
    assert(b1 == b2);
    assert(b1 != b3);
    print_test_result("BeadIndex comparison", b1 == b2 && b1 != b3);
    
    BeadIndex b_default;
    assert(static_cast<i32>(b_default) == -1);
    print_test_result("BeadIndex default construction", static_cast<i32>(b_default) == -1);
}



// Test 6: SimulationParameters
void test_simulation_parameters() {
    std::cout << "\n=== Testing SimulationParameters ===" << std::endl;
    
    // Default construction
    SimulationParameters params_default;
    print_test_result("Default construction compiles", true);
    
    // Valid construction
    SimulationParameters params{2.0, 100, 50, 3};
    assert(params.beta == 2.0);
    assert(params.num_beads == 100);
    assert(params.num_particles == 50);
    assert(params.num_dimensions == 3);
    assert(params.time_step == 2.0 / 100.0);
    print_test_result("SimulationParameters construction", params.time_step == 0.02);
    
    // Default dimension = 3
    SimulationParameters params_2d{2.0, 100, 50};
    assert(params_2d.num_dimensions == 3);
    print_test_result("Default num_dimensions = 3", params_2d.num_dimensions == 3);
    
    // Test that validation works (these should throw)
    bool caught_exception = false;
    try {
        SimulationParameters invalid{2.0, 0, 50};  // num_beads < 1
    } catch (const std::invalid_argument&) {
        caught_exception = true;
    }
    assert(caught_exception);
    print_test_result("Validation catches invalid num_beads", caught_exception);
    
    caught_exception = false;
    try {
        SimulationParameters invalid{2.0, 100, 0};  // num_particles < 1
    } catch (const std::invalid_argument&) {
        caught_exception = true;
    }
    assert(caught_exception);
    print_test_result("Validation catches invalid num_particles", caught_exception);
    
    caught_exception = false;
    try {
        SimulationParameters invalid{2.0, 100, 50, 4};  // num_dimensions > 3
    } catch (const std::invalid_argument&) {
        caught_exception = true;
    }
    assert(caught_exception);
    print_test_result("Validation catches invalid num_dimensions", caught_exception);
    
    caught_exception = false;
    try {
        SimulationParameters invalid{-1.0, 100, 50};  // beta <= 0
    } catch (const std::invalid_argument&) {
        caught_exception = true;
    }
    assert(caught_exception);
    print_test_result("Validation catches invalid beta", caught_exception);
    
    // Equality
    SimulationParameters params1{2.0, 100, 50, 3};
    SimulationParameters params2{2.0, 100, 50, 3};
    SimulationParameters params3{1.0, 100, 50, 3};
    assert(params1 == params2);
    assert(!(params1 == params3));
    print_test_result("SimulationParameters equality", params1 == params2 && !(params1 == params3));
}

// Test 7: MCStatistics
void test_mc_statistics() {
    std::cout << "\n=== Testing MCStatistics ===" << std::endl;
    
    MCStatistics stats;
    assert(stats.total_sweeps == 0);
    assert(stats.accepted_moves == 0);
    assert(stats.rejected_moves == 0);
    assert(stats.acceptance_rate == 0.0);
    print_test_result("MCStatistics default construction", true);
    
    // Record some moves
    stats.record_move(true);   // accepted
    stats.record_move(true);   // accepted
    stats.record_move(false);  // rejected
    
    assert(stats.total_sweeps == 3);
    assert(stats.accepted_moves == 2);
    assert(stats.rejected_moves == 1);
    assert(stats.acceptance_rate == 2.0 / 3.0);
    print_test_result("MCStatistics record_move and acceptance rate", 
                     stats.total_sweeps == 3 && stats.acceptance_rate == 2.0/3.0);
    
    // Reset
    stats.reset();
    assert(stats.total_sweeps == 0);
    assert(stats.accepted_moves == 0);
    assert(stats.rejected_moves == 0);
    assert(stats.acceptance_rate == 0.0);
    print_test_result("MCStatistics reset", stats.total_sweeps == 0);
    
    // Test with all rejections
    MCStatistics stats2;
    stats2.record_move(false);
    stats2.record_move(false);
    assert(stats2.acceptance_rate == 0.0);
    print_test_result("MCStatistics with all rejections", stats2.acceptance_rate == 0.0);
    
    // Test with all acceptances
    MCStatistics stats3;
    stats3.record_move(true);
    stats3.record_move(true);
    assert(stats3.acceptance_rate == 1.0);
    print_test_result("MCStatistics with all acceptances", stats3.acceptance_rate == 1.0);
}

// Test 8: Numeric and Floating Concepts
void test_concepts() {
    std::cout << "\n=== Testing Concepts ===" << std::endl;
    
    bool is_numeric_int = Numeric<int>;
    bool is_numeric_float = Numeric<float>;
    bool is_numeric_string = Numeric<std::string>;
    
    assert(is_numeric_int);
    assert(is_numeric_float);
    assert(!is_numeric_string);
    print_test_result("Numeric concept", is_numeric_int && is_numeric_float && !is_numeric_string);
    
    bool is_floating_int = Floating<int>;
    bool is_floating_float = Floating<float>;
    bool is_floating_double = Floating<double>;
    
    assert(!is_floating_int);
    assert(is_floating_float);
    assert(is_floating_double);
    print_test_result("Floating concept", !is_floating_int && is_floating_float && is_floating_double);
    
    bool is_real_float = Real<float>;
    bool is_real_double = Real<double>;
    bool is_real_long_double = Real<long double>;
    bool is_real_int = Real<int>;
    
    assert(is_real_float);
    assert(is_real_double);
    assert(is_real_long_double);
    assert(!is_real_int);
    print_test_result("Real concept", is_real_float && is_real_double && is_real_long_double && !is_real_int);
    
    // Vector3 concept
    bool is_vector = Vector3<Vec3d>;
    bool is_not_vector = Vector3<int>;
    
    assert(is_vector);
    assert(!is_not_vector);
    print_test_result("Vector3 concept", is_vector && !is_not_vector);
}

// Test 9: Formatter (if using C++20 format)
void test_formatter() {
    std::cout << "\n=== Testing Formatter ===" << std::endl;
    
#ifdef __cpp_lib_format
    try {
        Mass mass{42.5};
        std::string formatted = std::format("Mass: {}", mass);
        assert(formatted == "Mass: 42.5");
        print_test_result("Formatter for StrongType", formatted == "Mass: 42.5");
    } catch (const std::exception& e) {
        print_test_result("Formatter for StrongType", false);
        std::cout << "    Formatter error: " << e.what() << std::endl;
    }
#else
    std::cout << "    Skipping formatter test (C++20 <format> not available)" << std::endl;
    print_test_result("Formatter for StrongType (skipped)", true);
#endif
}

// Main test runner
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "  QSim Core Types Test Suite" << std::endl;
    std::cout << "========================================" << std::endl;
    
    try {
        test_strong_type_basic();
        test_type_safety();
        test_vec3_operations();
        test_indices();
        test_simulation_parameters();
        test_mc_statistics();
        test_concepts();
        test_formatter();
        
        std::cout << "\n========================================" << std::endl;
        std::cout << "  All tests completed!" << std::endl;
        std::cout << "========================================" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "\nERROR: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
