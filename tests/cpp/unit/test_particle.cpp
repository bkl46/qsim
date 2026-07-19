// tests/cpp/unit/test_particle.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "../../../include/qsim/core/particle.hpp"

using namespace qsim;
using namespace Catch::Matchers;

TEST_CASE("Particle construction and properties", "[particle]") {
    SECTION("Default constructor creates invalid particle") {
        Particle p;
        REQUIRE(p.mass() == 0.0);
        REQUIRE(p.spin() == 0.0);
        REQUIRE(p.statistics() == Particle::Statistics::Boltzmann);
        REQUIRE(p.charge() == 0.0);
        REQUIRE(p.is_valid() == false);
    }


    SECTION("Boson construction") {
        Particle p = Particle::boson(4.0026);  // helium-4
        REQUIRE_THAT(p.mass(), WithinRel(4.0026, 1e-6));
        REQUIRE(p.statistics() == Particle::Statistics::Boson);
        REQUIRE(p.spin() == 0.0);
        REQUIRE(p.charge() == 0.0);
        REQUIRE(p.is_valid() == true);
    }

    SECTION("Fermion construction") {
        Particle p = Particle::fermion(1.0);  // electron
        REQUIRE_THAT(p.mass(), WithinRel(1.0, 1e-6));
        REQUIRE(p.statistics() == Particle::Statistics::Fermion);
        REQUIRE(p.spin() == 0.5);
        REQUIRE(p.charge() == 0);
        REQUIRE(p.is_valid() == true);
    }


    SECTION("Electron construction") {
        Particle p = Particle::electron();  // electron
        REQUIRE_THAT(p.mass(), WithinRel(1.0, 1e-6));
        REQUIRE(p.statistics() == Particle::Statistics::Fermion);
        REQUIRE(p.spin() == 0.5);
        REQUIRE(p.charge() == -1.0);
        REQUIRE(p.is_valid() == true);
    }



    SECTION("Custom particle with charge") {
        Particle p(1.0, Particle::Statistics::Fermion, 0.5, -1.0);
        REQUIRE(p.is_valid() == true);
    }

    SECTION("Invalid particle detection") {
        Particle p(-1.0, Particle::Statistics::Fermion, 0.5, -1.0);
        REQUIRE(p.is_valid() == false);
        
        Particle p1(-1.0, Particle::Statistics::Fermion, 0, -1.0);
        REQUIRE(p.is_valid() == false);  // Fermion needs half spin
    }
}

TEST_CASE("Particle comparison operators", "[particle]") {


    SECTION("Equality comparison") {
        auto p1 = Particle::boson(4.0026);
        auto p2 = Particle::boson(4.0026);
        auto p3 = Particle::fermion(4.0026);
        
        REQUIRE(p1 == p2);
        REQUIRE_FALSE(p1 == p3);
    }

    SECTION("Inequality comparison") {
        auto p1 = Particle::boson(4.0026);
        auto p2 = Particle::boson(4.0026);
        auto p3 = Particle::fermion(4.0026);
        
        REQUIRE_FALSE(p1 != p2);
        REQUIRE(p1 != p3);
    }

    SECTION("Different masses are not equal") {
        auto p1 = Particle::boson(4.0026);
        auto p2 = Particle::boson(3.0);
        REQUIRE(p1 != p2);
    }
}

TEST_CASE("Particle statistics helpers", "[particle]") {
    SECTION("is_boson helper") {
        auto boson = Particle::boson(4.0);
        auto fermion = Particle::fermion(1.0);
        auto boltzmann = Particle(1.0, Particle::Statistics::Boltzmann);
        
        REQUIRE(boson.is_boson() == true);
        REQUIRE(fermion.is_boson() == false);
        REQUIRE(boltzmann.is_boson() == false);
    }

    SECTION("is_fermion helper") {
        auto boson = Particle::boson(4.0);
        auto fermion = Particle::fermion(1.0);
        auto boltzmann = Particle(1.0, Particle::Statistics::Boltzmann);
        
        REQUIRE(boson.is_fermion() == false);
        REQUIRE(fermion.is_fermion() == true);
        REQUIRE(boltzmann.is_fermion() == false);
    }

    SECTION("statistics to string") {
        auto boson = Particle::boson(4.0);
        auto fermion = Particle::fermion(1.0);
        auto boltzmann = Particle(1.0, Particle::Statistics::Boltzmann);
        
        REQUIRE(boson.statistics_string() == "Boson");
        REQUIRE(fermion.statistics_string() == "Fermion");
        REQUIRE(boltzmann.statistics_string() == "Boltzmann");
    }
}

TEST_CASE("Particle common species factories", "[particle]") {

    SECTION("Helium-4 boson") {
        auto he4 = Particle::helium4();
        REQUIRE_THAT(he4.mass(), WithinRel(7294.299, 1e-4));
        REQUIRE(he4.statistics() == Particle::Statistics::Boson);
        REQUIRE(he4.spin() == 0.0);
        REQUIRE(he4.charge() == 0.0);
    }

    SECTION("Electron fermion") {
        auto electron = Particle::electron();
        REQUIRE_THAT(electron.mass(), WithinRel(1.0, 1e-6));  // in a.u.
        REQUIRE(electron.statistics() == Particle::Statistics::Fermion);
        REQUIRE(electron.spin() == 0.5);
        REQUIRE_THAT(electron.charge(), WithinRel(-1.0, 1e-6));
    }

    SECTION("Proton fermion") {
        auto proton = Particle::proton();
        REQUIRE_THAT(proton.mass(), WithinRel(1836.15, 1e-3));
        REQUIRE(proton.statistics() == Particle::Statistics::Fermion);
        REQUIRE(proton.spin() == 0.5);
        REQUIRE_THAT(proton.charge(), WithinRel(1.0, 1e-6));
    }
}

TEST_CASE("Particle hashing for containers", "[particle]") {
    
    SECTION("Particle can be used in unordered containers") {
        auto p1 = Particle::boson(4.0026);
        auto p2 = Particle::fermion(1.0);
        
        std::unordered_map<Particle, double, ParticleHash> map;
        map[p1] = 4.0026;
        map[p2] = 1.0;
        
        REQUIRE(map.size() == 2);
        REQUIRE_THAT(map[p1], WithinRel(4.0026, 1e-6));
    }

    SECTION("Same particles hash to same value") {
        auto p1 = Particle::boson(4.0026);
        auto p2 = Particle::boson(4.0026);
        
        ParticleHash hasher;
        REQUIRE(hasher(p1) == hasher(p2));
    }

    SECTION("Different particles hash to different values") {
        auto p1 = Particle::boson(4.0026);
        auto p2 = Particle::fermion(4.0026);
        
        ParticleHash hasher;
        REQUIRE(hasher(p1) != hasher(p2));
    }

}

TEST_CASE("Particle copy and move semantics", "[particle]") {
    
    SECTION("Copy constructor") {
        auto original = Particle::boson(4.0026);
        auto copy = original;
        
        REQUIRE(copy.mass() == original.mass());
        REQUIRE(copy.statistics() == original.statistics());
        REQUIRE(copy.spin() == original.spin());
        REQUIRE(copy.charge() == original.charge());
    }

    SECTION("Copy assignment") {
        auto original = Particle::boson(4.0026);
        Particle copy;
        copy = original;
        
        REQUIRE(copy.mass() == original.mass());
        REQUIRE(copy.statistics() == original.statistics());
    }

    SECTION("Move constructor") {
        auto original = Particle::boson(4.0026);
        auto moved = std::move(original);
        
        REQUIRE_THAT(moved.mass(), WithinRel(4.0026, 1e-6));
        REQUIRE(moved.statistics() == Particle::Statistics::Boson);
    }
    
}

TEST_CASE("Particle physical properties", "[particle]") {
    
    SECTION("Compton wavelength scaling") {
        auto electron = Particle::electron();
        auto proton = Particle::proton();
        
        // Proton should have smaller de Broglie wavelength at same T
        // Compare reduced Compton wavelength: λ_c = hbar/(mc)
        double electron_lambda = 1.0 / electron.mass();
        double proton_lambda = 1.0 / proton.mass();
        
        REQUIRE(proton_lambda < electron_lambda);
    }

    SECTION("Mass comparison for sorting") {
        auto electron = Particle::electron();
        auto proton = Particle::proton();
        auto helium = Particle::helium4();


        
        std::vector<Particle> particles = {proton, electron, helium};
        std::sort(particles.begin(), particles.end(),
                  [](const Particle& a, const Particle& b) {
                      return a.mass() < b.mass();
                  });

        REQUIRE(particles[0].mass() == electron.mass());
        REQUIRE(particles[1].mass() == proton.mass());
        REQUIRE(particles[2].mass() == helium.mass());
    }
    
}
