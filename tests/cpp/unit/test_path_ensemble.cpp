// tests/unit/test_path_ensemble.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <random>

#include "qsim/path/path_ensemble.hpp"
#include "qsim/core/system.hpp"
#include "qsim/core/particle.hpp"
#include "qsim/potential/harmonic.hpp"

using namespace qsim;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

class PathEnsembleTestFixture {
public:
    PathEnsembleTestFixture() {
        system_.set_temperature(100.0);
        system_.set_num_beads(16);
        system_.set_box_size(5.0);
        system_.set_boundary(BoundaryType::Periodic);
        system_.set_dimensions(3);
        
        helium_ = Particle("He4", 4.0026, Statistics::Boson);
        system_.add_particles(helium_, 2);
        
        ensemble_ = std::make_unique<PathEnsemble>(system_);
    }
    
    System system_;
    Particle helium_;
    std::unique_ptr<PathEnsemble> ensemble_;
};

// ========================================================================
// Construction Tests
// ========================================================================

TEST_CASE("PathEnsemble construction", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    REQUIRE(ensemble.get_num_particles() == 2);
    REQUIRE(ensemble.get_num_beads() == 16);
    REQUIRE(ensemble.get_dimensions() == 3);
    REQUIRE(ensemble.data_size() == 2 * 16 * 3);
    REQUIRE_FALSE(ensemble.empty());
}

TEST_CASE("PathEnsemble copy and move", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    ensemble.set_bead_vector(0, 0, Vector3D{1.0, 2.0, 3.0});
    
    SECTION("Copy constructor") {
        PathEnsemble copy(ensemble);
        REQUIRE_THAT(copy.get_bead_vector(0, 0).x, WithinRel(1.0, 1e-12));
        REQUIRE_THAT(copy.get_bead_vector(0, 0).y, WithinRel(2.0, 1e-12));
        REQUIRE_THAT(copy.get_bead_vector(0, 0).z, WithinRel(3.0, 1e-12));
    }
    
    SECTION("Move constructor") {
        PathEnsemble moved(std::move(ensemble));
        REQUIRE_THAT(moved.get_bead_vector(0, 0).x, WithinRel(1.0, 1e-12));
        REQUIRE_THAT(moved.get_bead_vector(0, 0).y, WithinRel(2.0, 1e-12));
        REQUIRE_THAT(moved.get_bead_vector(0, 0).z, WithinRel(3.0, 1e-12));
    }
    
    SECTION("Copy assignment") {
        PathEnsemble assigned(fixture.system_);
        assigned = ensemble;
        REQUIRE_THAT(assigned.get_bead_vector(0, 0).x, WithinRel(1.0, 1e-12));
    }
    
    SECTION("Move assignment") {
        PathEnsemble assigned(fixture.system_);
        assigned = std::move(ensemble);
        REQUIRE_THAT(assigned.get_bead_vector(0, 0).x, WithinRel(1.0, 1e-12));
    }
}

// ========================================================================
// Accessor Tests
// ========================================================================

TEST_CASE("PathEnsemble bead access", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    SECTION("Set and get bead vector") {
        Vector3D pos{1.0, 2.0, 3.0};
        ensemble.set_bead_vector(0, 0, pos);
        
        Vector3D retrieved = ensemble.get_bead_vector(0, 0);
        REQUIRE_THAT(retrieved.x, WithinRel(1.0, 1e-12));
        REQUIRE_THAT(retrieved.y, WithinRel(2.0, 1e-12));
        REQUIRE_THAT(retrieved.z, WithinRel(3.0, 1e-12));
    }
    
    SECTION("Individual component access") {
        ensemble.set_bead(0, 0, 0, 1.0);
        ensemble.set_bead(0, 0, 1, 2.0);
        ensemble.set_bead(0, 0, 2, 3.0);
        
        REQUIRE_THAT(ensemble.get_bead(0, 0, 0), WithinRel(1.0, 1e-12));
        REQUIRE_THAT(ensemble.get_bead(0, 0, 1), WithinRel(2.0, 1e-12));
        REQUIRE_THAT(ensemble.get_bead(0, 0, 2), WithinRel(3.0, 1e-12));
        
        ensemble.set_bead(0, 0, 0, 5.0);
        REQUIRE_THAT(ensemble.get_bead(0, 0, 0), WithinRel(5.0, 1e-12));
    }
    
    SECTION("Invalid access throws") {
        REQUIRE_THROWS_AS(ensemble.get_bead(10, 0, 0), std::out_of_range);
        REQUIRE_THROWS_AS(ensemble.get_bead(0, 100, 0), std::out_of_range);
        REQUIRE_THROWS_AS(ensemble.get_bead(0, 0, 5), std::out_of_range);
    }
}

TEST_CASE("PathEnsemble set positions", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    std::vector<double> positions(96, 0.0);
    positions[0] = 1.0;
    positions[1] = 2.0;
    positions[2] = 3.0;
    
    ensemble.set_positions(positions);
    REQUIRE_THAT(ensemble.get_bead(0, 0, 0), WithinRel(1.0, 1e-12));
    REQUIRE_THAT(ensemble.get_bead(0, 0, 1), WithinRel(2.0, 1e-12));
    REQUIRE_THAT(ensemble.get_bead(0, 0, 2), WithinRel(3.0, 1e-12));
    
    SECTION("Invalid positions throw") {
        std::vector<double> wrong_size(10, 0.0);
        REQUIRE_THROWS_AS(ensemble.set_positions(wrong_size), std::invalid_argument);
    }
}

// ========================================================================
// Topology Tests
// ========================================================================

TEST_CASE("PathEnsemble topology management", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    SECTION("Set and get topology") {
        std::vector<ParticleID> topology = {0, 1};
        ensemble.set_topology(topology);
        REQUIRE(ensemble.get_topology().size() == 2);
        REQUIRE(ensemble.get_topology()[0] == 0);
        REQUIRE(ensemble.get_topology()[1] == 1);
    }
    
    SECTION("Swap paths") {
        ensemble.set_topology({0, 1});
        ensemble.swap_paths(0, 1);
        REQUIRE(ensemble.get_topology()[0] == 1);
        REQUIRE(ensemble.get_topology()[1] == 0);
    }
    
    SECTION("Invalid topology throws") {
        std::vector<ParticleID> wrong_size = {0, 1, 2};
        REQUIRE_THROWS_AS(ensemble.set_topology(wrong_size), std::invalid_argument);
    }
}

// ========================================================================
// Memory Management Tests
// ========================================================================

TEST_CASE("PathEnsemble memory management", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    SECTION("Clear") {
        ensemble.set_bead(0, 0, 0, 1.0);
        ensemble.clear();
        REQUIRE_THAT(ensemble.get_bead(0, 0, 0), WithinRel(0.0, 1e-12));
    }
    
    SECTION("Resize") {
        ensemble.resize(3, 32);
        REQUIRE(ensemble.get_num_particles() == 3);
        REQUIRE(ensemble.get_num_beads() == 32);
        REQUIRE(ensemble.data_size() == 3 * 32 * 3);
    }
}

// ========================================================================
// Physics Computation Tests
// ========================================================================

TEST_CASE("PathEnsemble kinetic action", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    SECTION("All beads at origin - zero action") {
        ensemble.clear();
        double action = ensemble.compute_kinetic_action();
        REQUIRE_THAT(action, WithinAbs(0.0, 1e-12));
    }
    
    SECTION("Random positions - positive action") {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        
        for (size_t i = 0; i < ensemble.data_size(); ++i) {
            ensemble.data()[i] = dist(gen);
        }
        
        double action = ensemble.compute_kinetic_action();
        REQUIRE(action > 0.0);
    }
}

TEST_CASE("PathEnsemble potential action", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& system = fixture.system_;
    auto& ensemble = *fixture.ensemble_;
    
    system.add_potential(std::make_unique<HarmonicPotential>(1.0, 0.0));
    
    SECTION("All beads at origin - zero action") {
        ensemble.clear();
        double action = ensemble.compute_potential_action();
        REQUIRE_THAT(action, WithinAbs(0.0, 1e-12));
    }
    
    SECTION("All beads at x=1.0") {
        for (size_t p = 0; p < ensemble.get_num_particles(); ++p) {
            for (uint32_t b = 0; b < ensemble.get_num_beads(); ++b) {
                ensemble.set_bead(p, b, 0, 1.0);
            }
        }
        
        double action = ensemble.compute_potential_action();
        // Each bead has energy 0.5, times 2 particles * 16 beads = 16.0
        REQUIRE_THAT(action, WithinRel(16.0, 1e-9));
    }
}

TEST_CASE("PathEnsemble action delta", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& system = fixture.system_;
    auto& ensemble = *fixture.ensemble_;
    
    system.add_potential(std::make_unique<HarmonicPotential>(1.0, 0.0));
    ensemble.clear();
    
    SECTION("Move bead from origin to (1,0,0)") {
        Vector3D new_pos{1.0, 0.0, 0.0};
        double delta = ensemble.compute_action_delta(0, 0, new_pos);
        
        // Should include kinetic + potential contributions
        REQUIRE(delta != 0.0);
        
        // Apply the update
        ensemble.update_bead(0, 0, new_pos);
        REQUIRE_THAT(ensemble.get_bead(0, 0, 0), WithinRel(1.0, 1e-12));
    }
}

TEST_CASE("PathEnsemble total action", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& system = fixture.system_;
    auto& ensemble = *fixture.ensemble_;
    
    system.add_potential(std::make_unique<HarmonicPotential>(1.0, 0.0));
    ensemble.clear();
    
    double kinetic = ensemble.compute_kinetic_action();
    double potential = ensemble.compute_potential_action();
    double total = ensemble.compute_total_action();
    
    REQUIRE_THAT(total, WithinRel(kinetic + potential, 1e-12));
}

// ========================================================================
// Modifier Tests
// ========================================================================

TEST_CASE("PathEnsemble update bead", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    ensemble.clear();
    Vector3D new_pos{1.0, 2.0, 3.0};
    ensemble.update_bead(0, 0, new_pos);
    
    REQUIRE_THAT(ensemble.get_bead(0, 0, 0), WithinRel(1.0, 1e-12));
    REQUIRE_THAT(ensemble.get_bead(0, 0, 1), WithinRel(2.0, 1e-12));
    REQUIRE_THAT(ensemble.get_bead(0, 0, 2), WithinRel(3.0, 1e-12));
}

TEST_CASE("PathEnsemble update multiple beads", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    ensemble.clear();
    std::vector<uint32_t> bead_indices = {0, 1, 2};
    std::vector<Vector3D> new_positions = {
        {1.0, 0.0, 0.0},
        {2.0, 0.0, 0.0},
        {3.0, 0.0, 0.0}
    };
    
    ensemble.update_beads(0, bead_indices, new_positions);
    
    REQUIRE_THAT(ensemble.get_bead(0, 0, 0), WithinRel(1.0, 1e-12));
    REQUIRE_THAT(ensemble.get_bead(0, 1, 0), WithinRel(2.0, 1e-12));
    REQUIRE_THAT(ensemble.get_bead(0, 2, 0), WithinRel(3.0, 1e-12));
}

TEST_CASE("PathEnsemble translate all", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    ensemble.clear();
    ensemble.set_bead_vector(0, 0, Vector3D{1.0, 0.0, 0.0});
    
    ensemble.translate_all(Vector3D{2.0, 3.0, 4.0});
    
    REQUIRE_THAT(ensemble.get_bead(0, 0, 0), WithinRel(3.0, 1e-12));
    REQUIRE_THAT(ensemble.get_bead(0, 0, 1), WithinRel(3.0, 1e-12));
    REQUIRE_THAT(ensemble.get_bead(0, 0, 2), WithinRel(4.0, 1e-12));
}

TEST_CASE("PathEnsemble periodic boundary", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    double box = 5.0;
    ensemble.clear();
    ensemble.set_bead(0, 0, 0, 6.0);
    
    ensemble.apply_periodic_boundary(box);
    
    REQUIRE_THAT(ensemble.get_bead(0, 0, 0), WithinRel(1.0, 1e-12));  // 6 - 5 = 1
}

// ========================================================================
// Estimator Tests
// ========================================================================

TEST_CASE("PathEnsemble center of mass", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    SECTION("All beads at (1,2,3)") {
        for (uint32_t b = 0; b < ensemble.get_num_beads(); ++b) {
            ensemble.set_bead_vector(0, b, Vector3D{1.0, 2.0, 3.0});
        }
        
        Vector3D com = ensemble.compute_path_center_of_mass(0);
        REQUIRE_THAT(com.x, WithinRel(1.0, 1e-12));
        REQUIRE_THAT(com.y, WithinRel(2.0, 1e-12));
        REQUIRE_THAT(com.z, WithinRel(3.0, 1e-12));
    }
}

TEST_CASE("PathEnsemble winding number", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    double box = 5.0;
    
    SECTION("No winding") {
        ensemble.clear();
        Vector3D winding = ensemble.compute_winding_number(0, box);
        REQUIRE_THAT(winding.x, WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(winding.y, WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(winding.z, WithinAbs(0.0, 1e-12));
    }
    
    SECTION("Path winds around box") {
        for (uint32_t b = 0; b < ensemble.get_num_beads(); ++b) {
            double frac = static_cast<double>(b) / ensemble.get_num_beads();
            ensemble.set_bead(0, b, 0, frac * box + box/2.0);
            ensemble.set_bead(0, b, 1, 0.0);
            ensemble.set_bead(0, b, 2, 0.0);
        }
        
        Vector3D winding = ensemble.compute_winding_number(0, box);
        REQUIRE(winding.x != 0.0);  // Should have winding
    }
}

TEST_CASE("PathEnsemble mean squared displacement", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    ensemble.clear();
    
    for (uint32_t b = 0; b < ensemble.get_num_beads(); ++b) {
        double pos = static_cast<double>(b);
        ensemble.set_bead(0, b, 0, pos);
    }
    
    double msd = ensemble.compute_mean_squared_displacement(0, 1);
    REQUIRE(msd > 0.0);
}

// ========================================================================
// Serialization Tests
// ========================================================================

TEST_CASE("PathEnsemble serialization", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& system = fixture.system_;
    auto& ensemble = *fixture.ensemble_;
    
    // Set some data
    ensemble.set_bead_vector(0, 0, Vector3D{1.0, 2.0, 3.0});
    ensemble.set_bead_vector(1, 0, Vector3D{4.0, 5.0, 6.0});
    
    const std::string filename = "test_ensemble.bin";
    
    SECTION("Save and load") {
        REQUIRE_NOTHROW(ensemble.save_to_binary(filename));
        
        PathEnsemble loaded = PathEnsemble::load_from_binary(filename, system);
        
        REQUIRE(loaded.get_num_particles() == ensemble.get_num_particles());
        REQUIRE(loaded.get_num_beads() == ensemble.get_num_beads());
        REQUIRE(loaded.get_dimensions() == ensemble.get_dimensions());
        
        auto v1 = loaded.get_bead_vector(0, 0);
        REQUIRE_THAT(v1.x, WithinRel(1.0, 1e-12));
        REQUIRE_THAT(v1.y, WithinRel(2.0, 1e-12));
        REQUIRE_THAT(v1.z, WithinRel(3.0, 1e-12));
        
        auto v2 = loaded.get_bead_vector(1, 0);
        REQUIRE_THAT(v2.x, WithinRel(4.0, 1e-12));
        REQUIRE_THAT(v2.y, WithinRel(5.0, 1e-12));
        REQUIRE_THAT(v2.z, WithinRel(6.0, 1e-12));
    }
    
    // Clean up
    std::remove(filename.c_str());
}

// ========================================================================
// Comparison Tests
// ========================================================================

TEST_CASE("PathEnsemble comparison", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& system = fixture.system_;
    auto& ensemble = *fixture.ensemble_;
    
    PathEnsemble same(system);
    same = ensemble;
    
    PathEnsemble different(system);
    different.set_bead(0, 0, 0, 999.0);
    
    SECTION("Equal ensembles") {
        REQUIRE(ensemble == same);
        REQUIRE_FALSE(ensemble != same);
    }
    
    SECTION("Different ensembles") {
        REQUIRE(ensemble != different);
        REQUIRE_FALSE(ensemble == different);
    }
}

// ========================================================================
// Validation Tests
// ========================================================================

TEST_CASE("PathEnsemble validation", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& ensemble = *fixture.ensemble_;
    
    SECTION("Valid ensemble") {
        ensemble.clear();
        REQUIRE(ensemble.is_valid());
    }
    
    SECTION("Invalid ensemble - NaN") {
        ensemble.set_bead(0, 0, 0, std::numeric_limits<double>::quiet_NaN());
        REQUIRE_FALSE(ensemble.is_valid());
    }
    
    SECTION("Invalid ensemble - topology mismatch") {
        // This would happen if someone manually corrupts the data
        // We'll test the size check
        std::vector<ParticleID> wrong_topo = {0, 1, 2};  // Should have 2
        REQUIRE_THROWS_AS(ensemble.set_topology(wrong_topo), std::invalid_argument);
    }
}

// ========================================================================
// Free Function Tests
// ========================================================================

TEST_CASE("PathEnsemble free functions", "[pathensemble]") {
    PathEnsembleTestFixture fixture;
    auto& system = fixture.system_;
    auto& ensemble = *fixture.ensemble_;
    
    SECTION("initialize_random") {
        initialize_random(ensemble, 5.0);
        
        bool all_zero = true;
        for (size_t i = 0; i < ensemble.data_size(); ++i) {
            if (ensemble.data()[i] != 0.0) {
                all_zero = false;
                break;
            }
        }
        REQUIRE_FALSE(all_zero);
    }
    
    SECTION("initialize_straight_paths") {
        ensemble.clear();
        Vector3D pos{1.0, 2.0, 3.0};
        ensemble.set_bead_vector(0, 0, pos);
        initialize_straight_paths(ensemble);
        
        // All beads should be at same position
        for (uint32_t b = 1; b < ensemble.get_num_beads(); ++b) {
            Vector3D bead_pos = ensemble.get_bead_vector(0, b);
            REQUIRE_THAT(bead_pos.x, WithinRel(pos.x, 1e-12));
            REQUIRE_THAT(bead_pos.y, WithinRel(pos.y, 1e-12));
            REQUIRE_THAT(bead_pos.z, WithinRel(pos.z, 1e-12));
        }
    }
    
    SECTION("bead_distance") {
        ensemble.clear();
        ensemble.set_bead_vector(0, 0, Vector3D{0.0, 0.0, 0.0});
        ensemble.set_bead_vector(0, 1, Vector3D{1.0, 0.0, 0.0});
        
        double dist = bead_distance(ensemble, 0, 0, 1, 5.0);
        REQUIRE_THAT(dist, WithinRel(1.0, 1e-12));
    }
    
    SECTION("kinetic_energy_between_beads") {
        double ke = kinetic_energy_between_beads(
            Vector3D{0.0, 0.0, 0.0},
            Vector3D{1.0, 0.0, 0.0},
            4.0, 0.1
        );
        REQUIRE(ke > 0.0);
    }
}
