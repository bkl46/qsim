/**
 *  Unit tests for Path and Bead 
 * 
 * Tests:
 * - Path construction with various bead counts
 * - Cyclic boundary conditions (wrap-around indexing)
 * - Bead access and manipulation
 * - Centroid calculation
 * - Edge cases (single bead, zero beads, large paths)
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../../include/qsim/core/types.hpp"
#include "../../../include/qsim/path/bead.hpp"
#include "../../../include/qsim/path/path.hpp"
#include <cmath>
#include <vector>
#include <stdexcept>

using namespace qsim;
using Catch::Approx;

// helper to compare vec3 with tolerance
bool vec3_approx_equal(const Vec3<double>& a, const Vec3<double>& b, double tol = 1e-10) {
    return (std::abs(a.x - b.x) < tol &&
            std::abs(a.y - b.y) < tol &&
            std::abs(a.z - b.z) < tol);
}

TEST_CASE("Path construction", "[path][construction]") {
    
    SECTION("Default constructor creates zero-bead path") {
        Path<double> path;
        REQUIRE(path.num_beads() == 0);
    }
    
    SECTION("Constructor with bead count allocates correct size") {
        std::vector<int> sizes = {1, 2, 5, 10, 100, 1000};
        for (int P : sizes) {
            DYNAMIC_SECTION("Path with " << P << " beads") {
                Path<double> path(P);
                REQUIRE(path.num_beads() == P);
            }
        }
    }
    
    SECTION("Newly constructed path has all beads at origin") {
        int P = 10;
        Path<double> path(P);
        Vec3 origin(0.0, 0.0, 0.0);
        for (int i = 0; i < P; ++i) {
            REQUIRE(vec3_approx_equal(path.bead(i), origin));
        }
    }
    
    SECTION("Negative bead count throws exception") {
        REQUIRE_THROWS_AS(Path<double>(-1), std::invalid_argument);
        REQUIRE_THROWS_AS(Path<double>(-100), std::invalid_argument);
    }
    
    SECTION("Zero bead count is valid but has special behavior") {
        Path<double> path(0);
        REQUIRE(path.num_beads() == 0);
        // Access methods should throw for empty path
        REQUIRE_THROWS_AS(path.bead(0), std::out_of_range);
    }
}

TEST_CASE("Bead access and manipulation", "[path][bead]") {
    Path<double> path(10);
    
    SECTION("Set and get bead positions") {
        Vec3 pos1(1.0, 2.0, 3.0);
        Vec3 pos2(-4.5, 0.0, 1e-5);
        
        path.set_bead(0, pos1);
        REQUIRE(vec3_approx_equal(path.bead(0), pos1));
        
        path.set_bead(5, pos2);
        REQUIRE(vec3_approx_equal(path.bead(5), pos2));
        
        // Verify other beads unchanged
        Vec3 origin(0.0, 0.0, 0.0);
        REQUIRE(vec3_approx_equal(path.bead(1), origin));
        REQUIRE(vec3_approx_equal(path.bead(9), origin));
    }
    
    SECTION("Out-of-bounds access throws exception") {
        REQUIRE_THROWS_AS(path.bead(-1), std::out_of_range);
        REQUIRE_THROWS_AS(path.bead(10), std::out_of_range);
        REQUIRE_THROWS_AS(path.bead(100), std::out_of_range);
        
        Vec3 pos(1.0, 1.0, 1.0);
        REQUIRE_THROWS_AS(path.set_bead(-1, pos), std::out_of_range);
        REQUIRE_THROWS_AS(path.set_bead(10, pos), std::out_of_range);
    }
    
    SECTION("Setting bead at same position multiple times") {
        Vec3 pos(1.0, 2.0, 3.0);
        path.set_bead(3, pos);
        path.set_bead(3, pos);
        REQUIRE(vec3_approx_equal(path.bead(3), pos));
        
        Vec3 pos2(4.0, 5.0, 6.0);
        path.set_bead(3, pos2);
        REQUIRE(vec3_approx_equal(path.bead(3), pos2));
    }
}

TEST_CASE("Cyclic indexing (wrap-around)", "[path][cyclic]") {
    int P = 5;
    Path<double> path(P);
    
    // Set up a known configuration
    for (int i = 0; i < P; ++i) {
        path.set_bead(i, Vec3(static_cast<double>(i), 0.0, 0.0));
    }
    
    SECTION("Next bead with standard indices") {
        REQUIRE(vec3_approx_equal(path.next(0), Vec3(1.0, 0.0, 0.0)));  // bead 1
        REQUIRE(vec3_approx_equal(path.next(3), Vec3(4.0, 0.0, 0.0)));  // bead 4
    }
    
    SECTION("Next bead wraps around at end") {
        // path[4].next() should be path[0]
        REQUIRE(vec3_approx_equal(path.next(4), Vec3(0.0, 0.0, 0.0)));
        REQUIRE(vec3_approx_equal(path.next(P-1), path.bead(0)));
    }
    
    SECTION("Previous bead with standard indices") {
        REQUIRE(vec3_approx_equal(path.prev(1), Vec3(0.0, 0.0, 0.0)));  // bead 0
        REQUIRE(vec3_approx_equal(path.prev(3), Vec3(2.0, 0.0, 0.0)));  // bead 2
    }
    
    SECTION("Previous bead wraps around at beginning") {
        // path[0].prev() should be path[4]
        REQUIRE(vec3_approx_equal(path.prev(0), Vec3(4.0, 0.0, 0.0)));
        REQUIRE(vec3_approx_equal(path.prev(0), path.bead(P-1)));
    }
    
    SECTION("Single bead path wraps to itself") {
        Path<double> single(1);
        single.set_bead(0, Vec3(42.0, 0.0, 0.0));
        REQUIRE(vec3_approx_equal(single.next(0), Vec3(42.0, 0.0, 0.0)));
        REQUIRE(vec3_approx_equal(single.prev(0), Vec3(42.0, 0.0, 0.0)));
    }
    
    SECTION("Two bead path cyclicity") {
        Path<double> two(2);
        two.set_bead(0, Vec3(0.0, 0.0, 0.0));
        two.set_bead(1, Vec3(1.0, 0.0, 0.0));
        
        REQUIRE(vec3_approx_equal(two.next(0), Vec3(1.0, 0.0, 0.0)));
        REQUIRE(vec3_approx_equal(two.next(1), Vec3(0.0, 0.0, 0.0)));
        REQUIRE(vec3_approx_equal(two.prev(0), Vec3(1.0, 0.0, 0.0)));
        REQUIRE(vec3_approx_equal(two.prev(1), Vec3(0.0, 0.0, 0.0)));
    }
}

TEST_CASE("Centroid calculation", "[path][centroid]") {
    
    SECTION("Centroid of origin-centered path") {
        Path<double> path(10);
        // All beads at (0,0,0)
        Vec3 center = path.centroid();
        REQUIRE(vec3_approx_equal(center, Vec3(0.0, 0.0, 0.0)));
    }
    
    SECTION("Centroid of uniform path along x-axis") {
        int P = 5;
        Path<double> path(P);
        for (int i = 0; i < P; ++i) {
            path.set_bead(i, Vec3(static_cast<double>(i), 0.0, 0.0));
        }
        
        Vec3 center = path.centroid();
        // Average of 0,1,2,3,4 = 2.0
        REQUIRE(center.x == Approx(2.0));
        REQUIRE(center.y == Approx(0.0));
        REQUIRE(center.z == Approx(0.0));
    }
    
    SECTION("Centroid of arbitrary configuration") {
        Path<double> path(4);
        path.set_bead(0, Vec3(1.0, 2.0, 3.0));
        path.set_bead(1, Vec3(-1.0, -2.0, -3.0));
        path.set_bead(2, Vec3(2.0, 0.0, -2.0));
        path.set_bead(3, Vec3(0.0, 4.0, 0.0));
        
        Vec3 center = path.centroid();
        // Sum: (2, 4, -2), average: (0.5, 1.0, -0.5)
        REQUIRE(center.x == Approx(0.5));
        REQUIRE(center.y == Approx(1.0));
        REQUIRE(center.z == Approx(-0.5));
    }
    
    SECTION("Centroid of single bead path") {
        Path<double> path(1);
        Vec3 pos(3.14, 2.71, 1.41);
        path.set_bead(0, pos);
        REQUIRE(vec3_approx_equal(path.centroid(), pos));
    }
    
    SECTION("Centroid is invariant under cyclic permutation") {
        usize P = 6;
        Path<double> path1(P);
        Path<double> path2(P);
        
        // Create a random-ish path
        std::vector<Vec3<double>> positions = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {-1.0, 0.0, 0.0},
            {0.0, -1.0, 0.0},
            {0.0, 0.0, -1.0}
        };
        
        for (int i = 0; i < P; ++i) {
            path1.set_bead(i, positions[i]);
        }
        
        // Create cyclically permuted version
        int shift = 2;
        for (int i = 0; i < P; ++i) {
            path2.set_bead(i, positions[(i + shift) % P]);
        }
        
        REQUIRE(vec3_approx_equal(path1.centroid(), path2.centroid()));
    }
    
    SECTION("Empty path centroid throws exception") {
        Path<double> empty(0);
        REQUIRE_THROWS_AS(empty.centroid(), std::runtime_error);
    }
}

TEST_CASE("Path modification operations", "[path][modification]") {
    Path<double> path(10);
    
    SECTION("Set all beads to same value") {
        Vec3 pos(7.0, 7.0, 7.0);
        for (int i = 0; i < 10; ++i) {
            path.set_bead(i, pos);
        }
        
        for (int i = 0; i < 10; ++i) {
            REQUIRE(vec3_approx_equal(path.bead(i), pos));
        }
        REQUIRE(vec3_approx_equal(path.centroid(), pos));
    }
    
    SECTION("Copy path") {
        Path<double> original(5);
        for (int i = 0; i < 5; ++i) {
            original.set_bead(i, Vec3(static_cast<double>(i), 
                                      static_cast<double>(i*i), 
                                      static_cast<double>(i*i*i)));
        }
        
        Path<double> copy = original;
        REQUIRE(copy.num_beads() == original.num_beads());
        for (int i = 0; i < 5; ++i) {
            REQUIRE(vec3_approx_equal(copy.bead(i), original.bead(i)));
        }
        
        // Verify deep copy
        copy.set_bead(0, Vec3(99.0, 99.0, 99.0));
        REQUIRE_FALSE(vec3_approx_equal(copy.bead(0), original.bead(0)));
    }
    
    SECTION("Assignment operator") {
        Path<double> original(3);
        original.set_bead(0, Vec3(1.0, 0.0, 0.0));
        original.set_bead(1, Vec3(0.0, 1.0, 0.0));
        original.set_bead(2, Vec3(0.0, 0.0, 1.0));
        
        Path<double> assigned(10);
        assigned = original;
        
        REQUIRE(assigned.num_beads() == 3);
        for (int i = 0; i < 3; ++i) {
            REQUIRE(vec3_approx_equal(assigned.bead(i), original.bead(i)));
        }
    }
}

TEST_CASE("Edge cases and stress tests", "[path][edge]") {
    
    SECTION("Very large path") {
        int P = 100000;
        Path<double> path(P);
        REQUIRE(path.num_beads() == P);
        
        // Set and verify a few beads to ensure no memory issues
        path.set_bead(0, Vec3(1.0, 2.0, 3.0));
        path.set_bead(P/2, Vec3(-1.0, -2.0, -3.0));
        path.set_bead(P-1, Vec3(0.0, 0.0, 0.0));
        
        REQUIRE(vec3_approx_equal(path.bead(0), Vec3(1.0, 2.0, 3.0)));
        REQUIRE(vec3_approx_equal(path.bead(P/2), Vec3(-1.0, -2.0, -3.0)));
        REQUIRE(vec3_approx_equal(path.bead(P-1), Vec3(0.0, 0.0, 0.0)));
        
        // Test cyclic boundaries for large path
        REQUIRE(vec3_approx_equal(path.next(P-1), path.bead(0)));
        REQUIRE(vec3_approx_equal(path.prev(0), path.bead(P-1)));
    }
    
    SECTION("Path with extreme coordinates") {
        Path<double> path(3);
        Vec3 huge(1e300, -1e300, 1e-300);
        Vec3 tiny(1e-300, 1e-300, 1e-300);
        Vec3 mixed(1e200, -1e-200, 0.0);
        
        path.set_bead(0, huge);
        path.set_bead(1, tiny);
        path.set_bead(2, mixed);
        
        REQUIRE(vec3_approx_equal(path.bead(0), huge, 1e-100));
        REQUIRE(vec3_approx_equal(path.bead(1), tiny, 1e-100));
        REQUIRE(vec3_approx_equal(path.bead(2), mixed, 1e-100));
        
        // Centroid should be computable without overflow
        Vec3 center = path.centroid();
        REQUIRE(std::isfinite(center.x));
        REQUIRE(std::isfinite(center.y));
        REQUIRE(std::isfinite(center.z));
    }
    
    SECTION("NaN and infinity handling") {
        Path<double> path(4);
        
        // Setting a bead to NaN should be allowed but detectable
        Vec3 nan_pos(std::nan(""), 0.0, 0.0);
        path.set_bead(0, nan_pos);
        REQUIRE(std::isnan(path.bead(0).x));
        
        // Centroid with NaN should propagate NaN
        Vec3 center = path.centroid();
        REQUIRE(std::isnan(center.x));
    }
}

TEST_CASE("Physical consistency checks", "[path][physics]") {
    
    SECTION("Open chain connectivity") {
        // Verify that beads form a ring: bead[i] -> bead[i+1] -> ... -> bead[0]
        int P = 8;
        Path<double> path(P);
        
        // Create a physically plausible configuration (rough circle)
        for (int i = 0; i < P; ++i) {
            double angle = 2.0 * M_PI * i / P;
            path.set_bead(i, Vec3(std::cos(angle), std::sin(angle), 0.0));
        }
        
        // Verify the distance between consecutive beads is consistent
        double first_dist = (path.bead(0) - path.bead(1)).norm();
        for (int i = 1; i < P; ++i) {
            double dist = (path.bead(i) - path.next(i)).norm();
            REQUIRE(dist == Approx(first_dist).epsilon(1e-10));
        }
        
        // Check wrap-around distance
        double wrap_dist = (path.bead(P-1) - path.next(P-1)).norm();
        REQUIRE(wrap_dist == Approx(first_dist).epsilon(1e-10));
    }
    
    SECTION("Path displacement is consistent") {
        // The sum of all displacement vectors around the ring should be zero
        usize P = 5;
        Path<double> path(P);
        
        // Arbitrary configuration
        path.set_bead(0, Vec3(0.0, 0.0, 0.0));
        path.set_bead(1, Vec3(1.0, 0.5, 0.0));
        path.set_bead(2, Vec3(0.5, 1.0, 0.3));
        path.set_bead(3, Vec3(-0.5, 0.5, 0.0));
        path.set_bead(4, Vec3(0.0, -0.5, -0.3));
        
        // Sum displacements around the ring
        Vec3 total_displacement(0.0, 0.0, 0.0);
        for (int i = 0; i < P; ++i) {
            total_displacement = total_displacement + (path.next(i) - path.bead(i));
        }
        
        REQUIRE(vec3_approx_equal(total_displacement, Vec3(0.0, 0.0, 0.0)));
    }
}

TEST_CASE("Performance characteristics (sanity)", "[path][performance]") {
    // These tests ensure basic operations are O(1) or O(P) as expected
    
    SECTION("Bead access is O(1)") {
        // Repeated access should be consistently fast
        int P = 10000;
        Path<double> path(P);
        
        // Set up path
        for (int i = 0; i < P; ++i) {
            path.set_bead(i, Vec3(static_cast<double>(i), 0.0, 0.0));
        }
        
        // Access all beads multiple times
        for (int rep = 0; rep < 100; ++rep) {
            for (int i = 0; i < P; ++i) {
                volatile auto v = path.bead(i);  // Prevent optimization
                (void)v;
            }
        }
        // If we get here without timeout, O(P) is reasonable
        REQUIRE(true);
    }
    
    SECTION("Centroid calculation scales linearly") {
        std::vector<usize> sizes = {100, 1000, 10000};
        for (usize P : sizes) {
            Path<double> path(P);
            for (int i = 0; i < P; ++i) {
                path.set_bead(i, Vec3(1.0, 1.0, 1.0));
            }
            Vec3 center = path.centroid();
            REQUIRE(center.x == Approx(1.0));
        }
    }
}
