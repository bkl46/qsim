// tests/cpp/unit/test_types.cpp
#include "../../../include/qsim/core/types.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <type_traits>

using namespace qsim;
using namespace Catch::Matchers;

// Test Fixture for common test data
struct TypesTestFixture {
    Mass m1{1.0};
    Mass m2{2.0};
    Mass m3{3.0};
    Vec3d v1{1.0, 2.0, 3.0};
    Vec3d v2{4.0, 5.0, 6.0};
    Vec3d v{3.0, 4.0, 0.0};
};

// Test 1: Strong Type Basic Operations
TEST_CASE_METHOD(TypesTestFixture, "Strong Type Basic Operations", "[strong_type]") {
    SECTION("Construction and value access") {
        Mass m_default;
        REQUIRE(m_default.value() == 0.0);
        REQUIRE(m1.value() == 1.0);
    }

    SECTION("Arithmetic operations") {
        // addition
        Mass m_sum = m1 + m2;
        REQUIRE(m_sum.value() == 3.0);
        
        // subtraction
        Mass m_diff = m3 - m1;
        REQUIRE(m_diff.value() == 2.0);
        
        // multiplication
        Mass m_mul = m1 * m2;
        REQUIRE(m_mul.value() == 2.0);
        
        // division
        Mass m_div = m3 / m1;
        REQUIRE(m_div.value() == 3.0);
    }

    SECTION("Scalar operations") {
        // Scalar multiplication (left)
        Mass m_scalar_left = 2.0 * m1;
        REQUIRE(m_scalar_left.value() == 2.0);
        
        // Scalar multiplication (right)
        Mass m_scalar_right = m1 * 2.0;
        REQUIRE(m_scalar_right.value() == 2.0);
    }

    SECTION("Implicit conversion") {
        // Implicit conversion 
        double raw_value = m1;
        REQUIRE(raw_value == 1.0);
    }
}

// Test 2: Type Safety
TEST_CASE("Type Safety with Strong Types", "[strong_type][type_safety]") {
    Mass mass{5.0};
    Length length{3.0};
    Time time{2.0};  

    SECTION("Same type operations compile") {
        Mass mass_result = mass + mass;
        REQUIRE(mass_result.value() == 10.0);
    }

    SECTION("Different strong types are distinct types") {
        bool same_type = std::is_same_v<Mass, Length>;
        REQUIRE_FALSE(same_type);
    }

    SECTION("StrongTypeConcept is satisfied") {
        bool satisfies_concept = StrongTypeConcept<Mass>;
        REQUIRE(satisfies_concept);
    }

    SECTION("Scalar multiplication with different types") {
        auto scaled_mass = mass * 2.0;
        REQUIRE(scaled_mass.value() == 10.0);
    }
}


// Test 3: Vec3 Operations
TEST_CASE_METHOD(TypesTestFixture, "Vec3 Operations", "[vec3]") {
    SECTION("Construction and access") {
        Vec3d v_scalar{2.0};  // explicit constructor
        REQUIRE(v1.x == 1.0);
        REQUIRE(v1.y == 2.0);
        REQUIRE(v1.z == 3.0);
        REQUIRE(v1[0] == 1.0);
        REQUIRE(v1[1] == 2.0);
        REQUIRE(v1[2] == 3.0);
    }

    SECTION("Arithmetic operations") {
        // Addition
        Vec3d v_sum = v1 + v2;
        REQUIRE(v_sum.x == 5.0);
        REQUIRE(v_sum.y == 7.0);
        REQUIRE(v_sum.z == 9.0);
        
        // Subtraction
        Vec3d v_diff = v2 - v1;
        REQUIRE(v_diff.x == 3.0);
        REQUIRE(v_diff.y == 3.0);
        REQUIRE(v_diff.z == 3.0);
    }

    SECTION("Scalar operations") {
        // Scalar multiplication
        Vec3d v_scaled = v1 * 2.0;
        REQUIRE(v_scaled.x == 2.0);
        REQUIRE(v_scaled.y == 4.0);
        REQUIRE(v_scaled.z == 6.0);
        
        // Scalar division
        Vec3d v_div = v2 / 2.0;
        REQUIRE(v_div.x == 2.0);
        REQUIRE(v_div.y == 2.5);
        REQUIRE(v_div.z == 3.0);
    }

    SECTION("In-place operators") {
        Vec3d v_plus_eq = v1;
        v_plus_eq += v2;
        REQUIRE(v_plus_eq.x == 5.0);
        REQUIRE(v_plus_eq.y == 7.0);
        REQUIRE(v_plus_eq.z == 9.0);
        
        Vec3d v_minus_eq = v2;
        v_minus_eq -= v1;
        REQUIRE(v_minus_eq.x == 3.0);
        REQUIRE(v_minus_eq.y == 3.0);
        REQUIRE(v_minus_eq.z == 3.0);
    }

    SECTION("Vector operations") {
        // Dot product
        double dot_result = v1.dot(v2);
        double expected_dot = 1.0*4.0 + 2.0*5.0 + 3.0*6.0;  // 4 + 10 + 18 = 32
        REQUIRE_THAT(dot_result, WithinAbs(expected_dot, 1e-12));
        
        // Cross product: v1 x v2 = (2*6 - 3*5, 3*4 - 1*6, 1*5 - 2*4) = (-3, 6, -3)
        Vec3d v_cross = v1.cross(v2);
        REQUIRE(v_cross.x == -3.0);
        REQUIRE(v_cross.y == 6.0);
        REQUIRE(v_cross.z == -3.0);
    }

    SECTION("Norm operations") {
        // Norm squared
        double norm_sq = v.norm_sq();
        REQUIRE(norm_sq == 25.0);
        
        // Norm
        double norm = v.norm();
        REQUIRE(norm == 5.0);
    }

    SECTION("Equality") {
        Vec3d v_equal{3.0, 4.0, 0.0};
        REQUIRE(v == v_equal);
        
        Vec3d v_not_equal{1.0, 2.0, 3.0};
        REQUIRE_FALSE(v == v_not_equal);
    }

    SECTION("Type aliases") {
        Vec3f v_float{1.0f, 2.0f, 3.0f};
        bool alias_works = std::is_same_v<decltype(v_float), Vec3<float>>;
        REQUIRE(alias_works);
    }
}

// Test 4: ParticleIndex and BeadIndex
TEST_CASE("ParticleIndex and BeadIndex", "[indices]") {
    SECTION("ParticleIndex operations") {
        ParticleIndex p1{5};
        ParticleIndex p2{5};
        ParticleIndex p3{10};
        
        // Conversions
        REQUIRE(static_cast<i32>(p1) == 5);
        REQUIRE(static_cast<usize>(p1) == 5);
        
        // Comparisons
        REQUIRE(p1 == p2);
        REQUIRE(p1 != p3);
        
        // Null checks
        REQUIRE(p1 != nullptr);
        
        ParticleIndex p_null;
        REQUIRE(p_null == nullptr);
    }

    SECTION("ParticleIndex increment") {
        ParticleIndex p_inc{0};
        ++p_inc;
        REQUIRE(static_cast<i32>(p_inc) == 1);
        
        ParticleIndex p_post{0};
        p_post++;
        REQUIRE(static_cast<i32>(p_post) == 1);
    }

    SECTION("BeadIndex operations") {
        BeadIndex b1{5};
        BeadIndex b2{5};
        BeadIndex b3{10};
        
        // Conversions
        REQUIRE(static_cast<i32>(b1) == 5);
        REQUIRE(static_cast<usize>(b1) == 5);
        
        // Comparisons
        REQUIRE(b1 == b2);
        REQUIRE(b1 != b3);
        
        // Default construction
        BeadIndex b_default;
        REQUIRE(static_cast<i32>(b_default) == -1);
    }
}

// Test 5: SimulationParameters
TEST_CASE("SimulationParameters", "[parameters]") {
    SECTION("Default construction") {
        SimulationParameters params_default;
        REQUIRE(params_default.beta == 0.0);  
        REQUIRE(params_default.num_beads == 0);
        REQUIRE(params_default.num_particles == 0);
        REQUIRE(params_default.num_dimensions == 3);
    }

    SECTION("Valid construction") {
        SimulationParameters params{2.0, 100, 50, 3};
        REQUIRE(params.beta == 2.0);
        REQUIRE(params.num_beads == 100);
        REQUIRE(params.num_particles == 50);
        REQUIRE(params.num_dimensions == 3);
        REQUIRE(params.time_step == 2.0 / 100.0);
    }

    SECTION("Default dimension = 3") {
        SimulationParameters params_2d{2.0, 100, 50};
        REQUIRE(params_2d.num_dimensions == 3);
    }

    SECTION("Validation catches invalid parameters") {
        // num_beads < 1
        REQUIRE_THROWS_AS((SimulationParameters{2.0, 0, 50}), std::invalid_argument);
        
        // num_particles < 1
        REQUIRE_THROWS_AS((SimulationParameters{2.0, 100, 0}), std::invalid_argument);
        
        // num_dimensions > 3
        REQUIRE_THROWS_AS((SimulationParameters{2.0, 100, 50, 4}), std::invalid_argument);
        
        // beta <= 0
        REQUIRE_THROWS_AS((SimulationParameters{-1.0, 100, 50}), std::invalid_argument);
    }

    SECTION("Equality comparisons") {
        SimulationParameters params1{2.0, 100, 50, 3};
        SimulationParameters params2{2.0, 100, 50, 3};
        SimulationParameters params3{1.0, 100, 50, 3};
        
        REQUIRE(params1 == params2);
        REQUIRE_FALSE(params1 == params3);
    }
}

// Test 6: MCStatistics
TEST_CASE("MCStatistics", "[statistics]") {
    SECTION("Default construction") {
        MCStatistics stats;
        REQUIRE(stats.total_sweeps == 0);
        REQUIRE(stats.accepted_moves == 0);
        REQUIRE(stats.rejected_moves == 0);
        REQUIRE(stats.acceptance_rate == 0.0);
    }

    SECTION("Recording moves") {
        MCStatistics stats;
        stats.record_move(true);   // accepted
        stats.record_move(true);   // accepted
        stats.record_move(false);  // rejected
        
        REQUIRE(stats.total_sweeps == 3);
        REQUIRE(stats.accepted_moves == 2);
        REQUIRE(stats.rejected_moves == 1);
        REQUIRE_THAT(stats.acceptance_rate, WithinAbs(2.0/3.0, 1e-12));
    }

    SECTION("Reset functionality") {
        MCStatistics stats;
        stats.record_move(true);
        stats.record_move(false);
        stats.reset();
        
        REQUIRE(stats.total_sweeps == 0);
        REQUIRE(stats.accepted_moves == 0);
        REQUIRE(stats.rejected_moves == 0);
        REQUIRE(stats.acceptance_rate == 0.0);
    }

    SECTION("Edge cases") {
        // All rejections
        MCStatistics stats2;
        stats2.record_move(false);
        stats2.record_move(false);
        REQUIRE(stats2.acceptance_rate == 0.0);
        
        // All acceptances
        MCStatistics stats3;
        stats3.record_move(true);
        stats3.record_move(true);
        REQUIRE(stats3.acceptance_rate == 1.0);
    }
}

// Test 7: Concepts
TEST_CASE("Type Concepts", "[concepts]") {
    SECTION("Numeric concept") {
        REQUIRE(Numeric<int>);
        REQUIRE(Numeric<float>);
        REQUIRE_FALSE(Numeric<std::string>);
    }

    SECTION("Floating concept") {
        REQUIRE_FALSE(Floating<int>);
        REQUIRE(Floating<float>);
        REQUIRE(Floating<double>);
    }

    SECTION("Real concept") {
        REQUIRE(Real<float>);
        REQUIRE(Real<double>);
        REQUIRE(Real<long double>);
        REQUIRE_FALSE(Real<int>);
    }

    SECTION("Vector3 concept") {
        REQUIRE(Vector3<Vec3d>);
        REQUIRE_FALSE(Vector3<int>);
    }
}

// Test 8: Formatter (Conditional)
TEST_CASE("Formatter for StrongType", "[formatter][!mayfail]") {
#ifdef __cpp_lib_format
    try {
        Mass mass{42.5};
        std::string formatted = std::format("Mass: {}", mass);
        REQUIRE(formatted == "Mass: 42.5");
    } catch (const std::exception& e) {
        FAIL("Formatter error: " << e.what());
    }
#else
    // Skip test if C++20 format not available
    // This test will be marked as skipped
    std::cout << "    Skipping formatter test (C++20 <format> not available)" << std::endl;
    SUCCEED("Test skipped - C++20 format not available");
#endif
}
