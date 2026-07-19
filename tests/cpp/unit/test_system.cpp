// tests/unit/test_system.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "qsim/core/system.hpp"
#include "qsim/core/particle.hpp"
#include "qsim/potential/harmonic.hpp"
#include "qsim/potential/lennard_jones.hpp"
#include "qsim/path/path_ensemble.hpp"

using namespace qsim;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

class SystemTestFixture {
public:
    SystemTestFixture() {
        helium_ = Particle("He4", 4.0026, Statistics::Boson, 0.0);
        electron_ = Particle("Electron", 0.511, Statistics::Fermion, -1.0);
        
        system_.set_temperature(300.0);
        system_.set_num_beads(64);
        system_.set_box_size(10.0);
        system_.set_boundary(BoundaryType::Periodic);
    }
    
    Particle helium_;
    Particle electron_;
    System system_;
};

// ========================================================================
// Construction and Basic Configuration Tests
// ========================================================================

TEST_CASE("System default construction", "[system]") {
    System default_system;
    REQUIRE(default_system.get_num_particles() == 0);
    REQUIRE(default_system.get_num_beads() == 64);
    REQUIRE_THAT(default_system.get_temperature(), WithinRel(300.0, 1e-9));
}

TEST_CASE("System parameterized construction", "[system]") {
    SECTION("Single parameter constructor") {
        System custom_system(100.0);
        REQUIRE_THAT(custom_system.get_temperature(), WithinRel(100.0, 1e-9));
        REQUIRE_THAT(custom_system.get_beta(), 
                     WithinRel(1.0 / (constants::k_B * 100.0), 1e-9));
    }
    
    SECTION("Two parameter constructor") {
        System custom_system2(200.0, 128);
        REQUIRE_THAT(custom_system2.get_temperature(), WithinRel(200.0, 1e-9));
        REQUIRE(custom_system2.get_num_beads() == 128);
    }
}

TEST_CASE("System particle management", "[system]") {
    SystemTestFixture fixture;
    auto& system = fixture.system_;
    auto& helium = fixture.helium_;
    auto& electron = fixture.electron_;
    
    SECTION("Add single particle") {
        ParticleID id1 = system.add_particle(helium);
        REQUIRE(id1 == 0);
        REQUIRE(system.get_num_particles() == 1);
    }
    
    SECTION("Add multiple particles") {
        auto ids = system.add_particles(electron, 5);
        REQUIRE(ids.size() == 5);
        REQUIRE(system.get_num_particles() == 5);
        REQUIRE(ids[0] == 0);
        REQUIRE(ids[4] == 4);
    }
    
    SECTION("Get particle by ID") {
        system.add_particle(helium);
        const Particle& p = system.get_particle(0);
        REQUIRE_THAT(p.mass, WithinRel(4.0026, 1e-6));
        REQUIRE(p.statistics == Statistics::Boson);
    }
    
    SECTION("Get particle with invalid ID throws") {
        REQUIRE_THROWS_AS(system.get_particle(0), std::out_of_range);
        REQUIRE_THROWS_AS(system.get_particle(100), std::out_of_range);
    }
    
    SECTION("Check identical particles") {
        system.add_particle(helium);
        system.add_particle(helium);
        system.add_particle(electron);
        
        REQUIRE(system.are_identical(0, 1));   // Both helium
        REQUIRE_FALSE(system.are_identical(0, 2)); // Helium vs electron
    }
}

// ========================================================================
// Potential Management Tests
// ========================================================================

TEST_CASE("System potential management", "[system]") {
    SystemTestFixture fixture;
    auto& system = fixture.system_;
    auto& helium = fixture.helium_;
    
    SECTION("Add harmonic potential") {
        auto harmonic = std::make_unique<HarmonicPotential>(1.0, 0.0);
        system.add_potential(std::move(harmonic));
        REQUIRE(system.get_potentials().size() == 1);
    }
    
    SECTION("Add Lennard-Jones with particle pairs") {
        auto lj = std::make_unique<LennardJonesPotential>(10.22, 2.556);
        system.add_particles(helium, 2);
        std::vector<std::pair<ParticleID, ParticleID>> pairs = {{0, 1}};
        system.add_potential(std::move(lj), pairs);
        REQUIRE(system.get_potentials().size() == 1);
    }
}

TEST_CASE("System potential evaluation", "[system]") {
    SystemTestFixture fixture;
    auto& system = fixture.system_;
    auto& helium = fixture.helium_;
    
    system.add_particle(helium);
    system.add_particle(helium);
    system.add_potential(std::make_unique<HarmonicPotential>(1.0, 0.0));
    
    SECTION("Both particles at origin") {
        std::vector<double> positions = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double energy = system.evaluate_potential(positions);
        REQUIRE_THAT(energy, WithinAbs(0.0, 1e-12));
    }
    
    SECTION("Both particles at x=1.0") {
        std::vector<double> positions = {1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
        double energy = system.evaluate_potential(positions);
        REQUIRE_THAT(energy, WithinRel(1.0, 1e-9));  // 0.5 * k * x^2 for each
    }
}

// ========================================================================
// Thermodynamic Parameter Tests
// ========================================================================

TEST_CASE("System temperature and beta", "[system]") {
    SystemTestFixture fixture;
    auto& system = fixture.system_;
    
    system.set_temperature(50.0);
    REQUIRE_THAT(system.get_temperature(), WithinRel(50.0, 1e-9));
    REQUIRE_THAT(system.get_beta(), WithinRel(1.0 / (constants::k_B * 50.0), 1e-9));
}

TEST_CASE("System imaginary time step", "[system]") {
    SystemTestFixture fixture;
    auto& system = fixture.system_;
    
    system.set_temperature(100.0);
    system.set_num_beads(100);
    
    double beta = 1.0 / (constants::k_B * 100.0);
    REQUIRE_THAT(system.get_tau(), WithinRel(beta / 100.0, 1e-9));
}

// ========================================================================
// Factory Method Tests
// ========================================================================

TEST_CASE("System create initial ensemble", "[system]") {
    SystemTestFixture fixture;
    auto& system = fixture.system_;
    auto& helium = fixture.helium_;
    
    system.add_particles(helium, 4);
    system.set_num_beads(32);
    system.set_box_size(5.0);
    
    PathEnsemble ensemble = system.create_initial_ensemble(InitializationType::Random);
    
    REQUIRE(ensemble.get_num_particles() == 4);
    REQUIRE(ensemble.get_num_beads() == 32);
    REQUIRE(ensemble.get_dimensions() == 3);
    REQUIRE_FALSE(ensemble.empty());
    REQUIRE(ensemble.data_size() == 4 * 32 * 3);
}

TEST_CASE("System create crystal ensemble", "[system]") {
    SystemTestFixture fixture;
    auto& system = fixture.system_;
    auto& helium = fixture.helium_;
    
    system.add_particles(helium, 8);
    system.set_box_size(4.0);
    
    PathEnsemble ensemble = system.create_crystal_ensemble(2.0);
    
    // In crystal initialization, all beads of a particle are at same position
    Vector3D pos0 = ensemble.get_bead_vector(0, 0);
    Vector3D pos1 = ensemble.get_bead_vector(0, 1);
    
    REQUIRE_THAT(pos0.x, WithinRel(pos1.x, 1e-12));
    REQUIRE_THAT(pos0.y, WithinRel(pos1.y, 1e-12));
    REQUIRE_THAT(pos0.z, WithinRel(pos1.z, 1e-12));
}

TEST_CASE("System create ensemble from positions", "[system]") {
    SystemTestFixture fixture;
    auto& system = fixture.system_;
    auto& helium = fixture.helium_;
    
    system.add_particles(helium, 2);
    system.set_num_beads(16);
    
    std::vector<double> positions(96, 0.0);
    positions[0] = 1.0;
    positions[1] = 2.0;
    positions[2] = 3.0;
    
    PathEnsemble ensemble = system.create_ensemble_from_positions(positions);
    
    REQUIRE_THAT(ensemble.get_bead(0, 0, 0), WithinRel(1.0, 1e-12));
    REQUIRE_THAT(ensemble.get_bead(0, 0, 1), WithinRel(2.0, 1e-12));
    REQUIRE_THAT(ensemble.get_bead(0, 0, 2), WithinRel(3.0, 1e-12));
}

// ========================================================================
// Serialization Tests
// ========================================================================

TEST_CASE("System JSON serialization", "[system]") {
    SystemTestFixture fixture;
    auto& system = fixture.system_;
    auto& helium = fixture.helium_;
    
    system.add_particles(helium, 3);
    system.add_potential(std::make_unique<HarmonicPotential>(0.5, 0.0));
    system.set_temperature(77.0);
    system.set_box_size(15.0);
    system.set_dimensions(3);
    system.set_boundary(BoundaryType::Periodic);
    
    const std::string filename = "test_system.json";
    
    REQUIRE_NOTHROW(system.save_to_json(filename));
    
    System loaded_system = System::load_from_json(filename);
    
    REQUIRE(loaded_system.get_num_particles() == 3);
    REQUIRE_THAT(loaded_system.get_temperature(), WithinRel(77.0, 1e-9));
    REQUIRE_THAT(loaded_system.get_box_size(), WithinRel(15.0, 1e-9));
    REQUIRE(loaded_system.get_dimensions() == 3);
    REQUIRE(loaded_system.get_boundary() == BoundaryType::Periodic);
    
    // Clean up
    std::remove(filename.c_str());
}

// ========================================================================
// Validation Tests
// ========================================================================

TEST_CASE("System validation", "[system]") {
    SystemTestFixture fixture;
    auto& system = fixture.system_;
    auto& helium = fixture.helium_;
    
    SECTION("Empty system is invalid") {
        REQUIRE_FALSE(system.is_valid());
        REQUIRE_FALSE(system.get_validation_errors().empty());
    }
    
    SECTION("System with particles is valid") {
        system.add_particles(helium, 2);
        REQUIRE(system.is_valid());
    }
    
    SECTION("Invalid temperature makes system invalid") {
        system.add_particles(helium, 2);
        system.set_temperature(0.0);
        REQUIRE_FALSE(system.is_valid());
    }
}

// ========================================================================
// Free Function Tests
// ========================================================================

TEST_CASE("System create single particle system", "[system]") {
    Particle he("He4", 4.0, Statistics::Boson);
    auto sys = create_single_particle_system<System>(he, 50.0);
    
    REQUIRE(sys.get_num_particles() == 1);
    REQUIRE_THAT(sys.get_temperature(), WithinRel(50.0, 1e-9));
    
    const auto& p = sys.get_particle(0);
    REQUIRE_THAT(p.mass, WithinRel(4.0, 1e-6));
    REQUIRE(p.statistics == Statistics::Boson);
}

TEST_CASE("System create fluid system", "[system]") {
    Particle he("He4", 4.0, Statistics::Boson);
    auto sys = create_fluid_system<System>(he, 64, 0.14, 2.17);
    
    REQUIRE(sys.get_num_particles() == 64);
    REQUIRE_THAT(sys.get_temperature(), WithinRel(2.17, 1e-9));
    REQUIRE(sys.get_box_size() > 0.0);  // Box size computed from density
}
