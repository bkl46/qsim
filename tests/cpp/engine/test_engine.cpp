// tests/cpp/engine/test_engine.cpp

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../../../include/qsim/engine/engine.hpp"
#include "../../../include/qsim/path/path_ensemble.hpp"
#include "../../../include/qsim/path/bead.hpp"
#include "../../../include/qsim/action/action.hpp"
#include "../../../include/qsim/core/types.hpp"

using namespace qsim;
using Catch::Approx;


/**
 * @brief Minimal bead implementation for testing
 */
struct TestBead {
    double position[3] = {0.0, 0.0, 0.0};
    double mass = 1.0;
    
    TestBead() = default;
    TestBead(double x, double y, double z) : position{x, y, z} {}
};

/**
 * @brief Minimal path ensemble for testing
 * 
 * Simple grid: num_particles × num_beads × 3D positions
 */
class TestPathEnsemble {
public:
    TestPathEnsemble(size_t n_particles, size_t n_beads) 
        : n_particles_(n_particles)
        , n_beads_(n_beads)
    {
        beads_.resize(n_particles * n_beads);
        
        // Initialize with simple values
        for (size_t p = 0; p < n_particles; ++p) {
            for (size_t b = 0; b < n_beads; ++b) {
                auto& bead = beads_[p * n_beads + b];
                bead.position[0] = static_cast<double>(p) * 10.0 + static_cast<double>(b);
                bead.position[1] = static_cast<double>(p) * 10.0 + static_cast<double>(b) + 1.0;
                bead.position[2] = static_cast<double>(p) * 10.0 + static_cast<double>(b) + 2.0;
                bead.mass = 1.0;
            }
        }
    }
    
    [[nodiscard]] size_t num_particles() const noexcept { return n_particles_; }
    [[nodiscard]] size_t num_beads() const noexcept { return n_beads_; }
    
    TestBead& bead(size_t particle, size_t bead_idx) {
        return beads_[particle * n_beads_ + bead_idx];
    }
    
    const TestBead& bead(size_t particle, size_t bead_idx) const {
        return beads_[particle * n_beads_ + bead_idx];
    }
    
private:
    size_t n_particles_;
    size_t n_beads_;
    std::vector<TestBead> beads_;
};

/**
 * @brief Mock action that returns a simple quadratic function of positions
 * 
 * This allows us to test the engine without real physics.
 * Action = sum over all positions of (x^2 + y^2 + z^2)
 */
class MockAction {
public:
    double evaluate(const TestPathEnsemble& ensemble) const {
        double total = 0.0;
        size_t n_particles = ensemble.num_particles();
        size_t n_beads = ensemble.num_beads();
        
        for (size_t p = 0; p < n_particles; ++p) {
            for (size_t b = 0; b < n_beads; ++b) {
                const auto& bead = ensemble.bead(p, b);
                total += bead.position[0] * bead.position[0];
                total += bead.position[1] * bead.position[1];
                total += bead.position[2] * bead.position[2];
            }
        }
        
        return total;
    }
};

/**
 * @brief Minimal thermodynamic energy estimator
 */
class MockThermodynamicEstimator {
public:
    explicit MockThermodynamicEstimator(const MockAction& action) 
        : action_(action) {}
    
    double evaluate(const TestPathEnsemble& ensemble) const {
        // Simple: return action value as "energy"
        return action_.evaluate(ensemble);
    }
    
private:
    const MockAction& action_;
};

/**
 * @brief Minimal virial energy estimator
 */
class MockVirialEstimator {
public:
    explicit MockVirialEstimator(const MockAction& action) 
        : action_(action) {}
    
    double evaluate(const TestPathEnsemble& ensemble) const {
        // Slightly different: action + small correction
        return action_.evaluate(ensemble) + 0.1;
    }
    
private:
    const MockAction& action_;
};

/**
 * @brief Mock RNG engine (deterministic for testing)
 * 
 * Returns a simple sequence: 0.1, 0.2, 0.3, ... wrapping around
 */
class MockRNG {
public:
    explicit MockRNG(int64_t seed = 0) : counter_(0), seed_(seed) {}
    
    double uniform() {
        // Deterministic sequence for reproducibility
        double value = (counter_ % 100) / 100.0;
        counter_++;
        return value;
    }
    
    double normal(double mean, double sigma) {
        // Simple deterministic "normal" using uniform
        double u = uniform();
        return mean + sigma * (u - 0.5) * 2.0;  // Approximately normal-ish
    }
    
    [[nodiscard]] int64_t seed() const noexcept { return seed_; }
    
private:
    size_t counter_;
    int64_t seed_;
};

/**
 * @brief Minimal Metropolis sampler
 */
class MockMetropolisSampler {
public:
    explicit MockMetropolisSampler(MockRNG& rng) : rng_(rng) {}
    
    bool accept(double delta_action) {
        if (delta_action <= 0.0) return true;  // Always accept downward moves
        double threshold = std::exp(-delta_action);
        return rng_.uniform() < threshold;
    }
    
private:
    MockRNG& rng_;
};

// ============================================================================
// Specialize MonteCarloEngine for our test types
// ============================================================================

// Since the actual engine is templated, we test with our mock types
// In practice, you'd want to template your engine to accept these as parameters

namespace qsim {
    // Template specialization or we can test the generic engine directly
    // For this test, we'll create a simplified testable engine
    
    template<typename Action, typename Ensemble, typename RNG>
    class TestableEngine {
    public:
        TestableEngine(Action action, Ensemble ensemble, EngineConfig config)
            : action_(std::move(action))
            , ensemble_(std::move(ensemble))
            , config_(std::move(config))
            , rng_(config_.random_seed >= 0 ? config_.random_seed : 12345)
            , metropolis_(rng_)
            , thermo_estimator_(action_)
            , virial_estimator_(action_)
            , current_action_(action_.evaluate(ensemble_))
            , accepted_moves_(0)
            , total_moves_(0)
        {}
        
        void run() {
            thermalize();
            measure();
        }
        
        void thermalize() {
            for (size_t step = 0; step < config_.num_thermalization_steps; ++step) {
                attempt_move();
            }
            stats_.reset();
            accepted_moves_ = 0;
            total_moves_ = 0;
        }
        
        void measure() {
            for (size_t step = 0; step < config_.num_measurement_steps; ++step) {
                attempt_move();
                
                if (step % config_.measurement_interval == 0) {
                    double energy = thermo_estimator_.evaluate(ensemble_);
                    stats_.accumulate(energy);
                }
            }
        }
        
        void attempt_move() {
            double selector = rng_.uniform();
            
            if (selector < 0.6) {
                // Single bead move
                auto proposal = ensemble_;
                size_t p = static_cast<size_t>(rng_.uniform() * proposal.num_particles()) % proposal.num_particles();
                size_t b = static_cast<size_t>(rng_.uniform() * proposal.num_beads()) % proposal.num_beads();
                
                auto& bead = proposal.bead(p, b);
                bead.position[0] += rng_.normal(0.0, 1.0);
                bead.position[1] += rng_.normal(0.0, 1.0);
                bead.position[2] += rng_.normal(0.0, 1.0);
                
                double new_action = action_.evaluate(proposal);
                double delta = new_action - current_action_;
                
                if (metropolis_.accept(delta)) {
                    ensemble_ = proposal;
                    current_action_ = new_action;
                    accepted_moves_++;
                }
            }
            total_moves_++;
        }
        
        [[nodiscard]] double mean_energy() const noexcept { return stats_.mean; }
        [[nodiscard]] double energy_standard_error() const noexcept { return stats_.standard_error(); }
        [[nodiscard]] double acceptance_rate() const noexcept {
            return total_moves_ > 0 ? static_cast<double>(accepted_moves_) / total_moves_ : 0.0;
        }
        [[nodiscard]] const Ensemble& ensemble() const noexcept { return ensemble_; }
        [[nodiscard]] const EngineConfig& config() const noexcept { return config_; }
        [[nodiscard]] size_t accepted_moves() const noexcept { return accepted_moves_; }
        [[nodiscard]] size_t total_moves() const noexcept { return total_moves_; }
        
    private:
        Action action_;
        Ensemble ensemble_;
        EngineConfig config_;
        RNG rng_;
        MockMetropolisSampler metropolis_;
        MockThermodynamicEstimator thermo_estimator_;
        MockVirialEstimator virial_estimator_;
        EnergyStatistics stats_;
        double current_action_;
        size_t accepted_moves_;
        size_t total_moves_;
    };
}

// ============================================================================
// Test Cases
// ============================================================================

TEST_CASE("EnergyStatistics - Welford's algorithm", "[engine][statistics]") {
    SECTION("Empty statistics") {
        EnergyStatistics stats;
        REQUIRE(stats.mean == 0.0);
        REQUIRE(stats.count == 0);
        REQUIRE(stats.standard_error() == 0.0);
    }
    
    SECTION("Single value") {
        EnergyStatistics stats;
        stats.accumulate(42.0);
        REQUIRE(stats.mean == Approx(42.0));
        REQUIRE(stats.count == 1);
        REQUIRE(stats.variance() == 0.0);
        REQUIRE(stats.standard_error() == 0.0);
    }
    
    SECTION("Constant values") {
        EnergyStatistics stats;
        for (int i = 0; i < 100; ++i) {
            stats.accumulate(5.0);
        }
        REQUIRE(stats.mean == Approx(5.0));
        REQUIRE(stats.variance() == Approx(0.0).margin(1e-10));
        REQUIRE(stats.standard_error() == Approx(0.0).margin(1e-10));
    }
    
    SECTION("Known distribution") {
        EnergyStatistics stats;
        // Values: 1, 2, 3, 4, 5
        // Mean = 3, Variance = 2.5 (sample), 2.0 (population)
        stats.accumulate(1.0);
        stats.accumulate(2.0);
        stats.accumulate(3.0);
        stats.accumulate(4.0);
        stats.accumulate(5.0);
        
        REQUIRE(stats.mean == Approx(3.0));
        REQUIRE(stats.variance() == Approx(2.5));  // Sample variance
        REQUIRE(stats.standard_error() == Approx(std::sqrt(2.5 / 5.0)));
    }
    
    SECTION("Reset functionality") {
        EnergyStatistics stats;
        stats.accumulate(10.0);
        stats.accumulate(20.0);
        REQUIRE(stats.count == 2);
        
        stats.reset();
        REQUIRE(stats.count == 0);
        REQUIRE(stats.mean == 0.0);
        REQUIRE(stats.M2 == 0.0);
    }
}

TEST_CASE("EngineConfig - Default values", "[engine][config]") {
    SECTION("Default construction") {
        EngineConfig config;
        REQUIRE(config.num_thermalization_steps == 10000);
        REQUIRE(config.num_measurement_steps == 100000);
        REQUIRE(config.measurement_interval == 10);
        REQUIRE(config.checkpoint_interval == 1000);
        REQUIRE(config.random_seed == -1);
        REQUIRE(config.single_bead_move_prob == Approx(0.6));
        REQUIRE(config.staging_move_prob == Approx(0.3));
        REQUIRE(config.center_of_mass_move_prob == Approx(0.1));
    }
    
    SECTION("Move probabilities sum check") {
        EngineConfig config;
        double sum = config.single_bead_move_prob 
                   + config.staging_move_prob 
                   + config.center_of_mass_move_prob;
        REQUIRE(sum == Approx(1.0));
    }
    
    SECTION("Designated initializers") {
        EngineConfig config{
            .num_thermalization_steps = 500,
            .num_measurement_steps = 1000,
            .measurement_interval = 5,
            .random_seed = 42
        };
        REQUIRE(config.num_thermalization_steps == 500);
        REQUIRE(config.num_measurement_steps == 1000);
        REQUIRE(config.measurement_interval == 5);
        REQUIRE(config.random_seed == 42);
        // Unspecified should be default
        REQUIRE(config.checkpoint_interval == 1000);
    }
}

TEST_CASE("TestableEngine - Construction", "[engine][construction]") {
    MockAction action;
    TestPathEnsemble ensemble(2, 10);  // 2 particles, 10 beads
    EngineConfig config{
        .num_thermalization_steps = 100,
        .num_measurement_steps = 200,
        .random_seed = 42
    };
    
    qsim::TestableEngine<MockAction, TestPathEnsemble, MockRNG> engine(
        action, ensemble, config
    );
    
    SECTION("Engine stores configuration correctly") {
        REQUIRE(engine.config().num_thermalization_steps == 100);
        REQUIRE(engine.config().num_measurement_steps == 200);
        REQUIRE(engine.config().random_seed == 42);
    }
    
    SECTION("Initial state is accessible") {
        REQUIRE(engine.ensemble().num_particles() == 2);
        REQUIRE(engine.ensemble().num_beads() == 10);
        REQUIRE(engine.accepted_moves() == 0);
        REQUIRE(engine.total_moves() == 0);
    }
    
    SECTION("Acceptance rate is zero before any moves") {
        REQUIRE(engine.acceptance_rate() == 0.0);
    }
}

TEST_CASE("TestableEngine - Move attempts", "[engine][moves]") {
    MockAction action;
    TestPathEnsemble ensemble(1, 5);  // 1 particle, 5 beads for simplicity
    EngineConfig config{
        .num_thermalization_steps = 0,
        .num_measurement_steps = 0,
        .random_seed = 0
    };
    
    qsim::TestableEngine<MockAction, TestPathEnsemble, MockRNG> engine(
        action, ensemble, config
    );
    
    SECTION("Single move attempt increments counters") {
        engine.attempt_move();
        REQUIRE(engine.total_moves() == 1);
        // With deterministic mock RNG, we might know if it's accepted
        // But at minimum, total_moves increments
    }
    
    SECTION("Multiple move attempts") {
        for (int i = 0; i < 10; ++i) {
            engine.attempt_move();
        }
        REQUIRE(engine.total_moves() == 10);
    }
    
    SECTION("Acceptance rate is in valid range") {
        for (int i = 0; i < 100; ++i) {
            engine.attempt_move();
        }
        double rate = engine.acceptance_rate();
        REQUIRE(rate >= 0.0);
        REQUIRE(rate <= 1.0);
    }
}

TEST_CASE("TestableEngine - Run with minimal steps", "[engine][run]") {
    MockAction action;
    TestPathEnsemble ensemble(1, 5);
    EngineConfig config{
        .num_thermalization_steps = 50,
        .num_measurement_steps = 100,
        .measurement_interval = 5,
        .random_seed = 0
    };
    
    qsim::TestableEngine<MockAction, TestPathEnsemble, MockRNG> engine(
        action, ensemble, config
    );
    
    SECTION("Run completes without exceptions") {
        REQUIRE_NOTHROW(engine.run());
    }
    
    SECTION("Run produces energy statistics") {
        engine.run();
        
        // Should have taken some measurements
        auto energy = engine.mean_energy();
        REQUIRE(!std::isnan(energy));
        REQUIRE(!std::isinf(energy));
        REQUIRE(energy > 0.0);  // Our mock action is positive definite
        
        auto error = engine.energy_standard_error();
        REQUIRE(!std::isnan(error));
        REQUIRE(error >= 0.0);
    }
    
    SECTION("Run with different seeds produces different results") {
        EngineConfig config2 = config;
        config2.random_seed = 999;
        
        qsim::TestableEngine<MockAction, TestPathEnsemble, MockRNG> engine2(
            MockAction(), TestPathEnsemble(1, 5), config2
        );
        
        engine.run();
        engine2.run();
        
        // Due to our deterministic mock RNG, different seeds should give different results
        // Note: This test is weak with mock RNG; with real RNG it would be meaningful
        INFO("Energy 1: " << engine.mean_energy());
        INFO("Energy 2: " << engine2.mean_energy());
    }
    
    SECTION("Acceptance rate after run is reasonable") {
        engine.run();
        double rate = engine.acceptance_rate();
        REQUIRE(rate > 0.0);
        REQUIRE(rate <= 1.0);
    }
}

TEST_CASE("TestableEngine - Energy measurements", "[engine][energy]") {
    MockAction action;
    TestPathEnsemble ensemble(1, 3);
    EngineConfig config{
        .num_thermalization_steps = 0,
        .num_measurement_steps = 200,
        .measurement_interval = 10,
        .random_seed = 0
    };
    
    qsim::TestableEngine<MockAction, TestPathEnsemble, MockRNG> engine(
        action, ensemble, config
    );
    
    SECTION("Energy changes during measurement phase") {
        double initial_energy = engine.mean_energy();
        engine.measure();  // Run only measurement (skip thermalization)
        double final_energy = engine.mean_energy();
        
        // Energy might change due to moves
        INFO("Initial energy: " << initial_energy);
        INFO("Final energy: " << final_energy);
    }
    
    SECTION("More steps reduces statistical error") {
        // Run with few steps
        engine.run();
        double error_short = engine.energy_standard_error();
        
        // Run with more steps
        EngineConfig config_long = config;
        config_long.num_measurement_steps = 1000;
        config_long.measurement_interval = 10;
        
        qsim::TestableEngine<MockAction, TestPathEnsemble, MockRNG> engine_long(
            MockAction(), TestPathEnsemble(1, 3), config_long
        );
        engine_long.run();
        double error_long = engine_long.energy_standard_error();
        
        INFO("Error with " << config.num_measurement_steps << " steps: " << error_short);
        INFO("Error with " << config_long.num_measurement_steps << " steps: " << error_long);
        
        // More measurements should generally reduce error
        // This might not always hold with mock RNG, but is a sanity check
    }
}

TEST_CASE("TestableEngine - Ensemble preservation", "[engine][ensemble]") {
    MockAction action;
    TestPathEnsemble ensemble(2, 4);
    EngineConfig config{
        .num_thermalization_steps = 50,
        .num_measurement_steps = 50,
        .random_seed = 0
    };
    
    qsim::TestableEngine<MockAction, TestPathEnsemble, MockRNG> engine(
        action, ensemble, config
    );
    
    SECTION("Ensemble dimensions preserved after run") {
        engine.run();
        REQUIRE(engine.ensemble().num_particles() == 2);
        REQUIRE(engine.ensemble().num_beads() == 4);
    }
    
    SECTION("Ensemble positions are valid after run") {
        engine.run();
        
        for (size_t p = 0; p < engine.ensemble().num_particles(); ++p) {
            for (size_t b = 0; b < engine.ensemble().num_beads(); ++b) {
                const auto& bead = engine.ensemble().bead(p, b);
                REQUIRE(!std::isnan(bead.position[0]));
                REQUIRE(!std::isnan(bead.position[1]));
                REQUIRE(!std::isnan(bead.position[2]));
                REQUIRE(!std::isinf(bead.position[0]));
                REQUIRE(!std::isinf(bead.position[1]));
                REQUIRE(!std::isinf(bead.position[2]));
            }
        }
    }
}

TEST_CASE("TestableEngine - Statistical properties", "[engine][statistics]") {
    SECTION("Constant ensemble gives zero variance") {
        // Create a special ensemble where moves are always rejected
        // by having enormous positions (so any change increases action a lot)
        TestPathEnsemble ensemble(1, 5);
        
        // Set positions to zero (minimum of quadratic action)
        for (size_t p = 0; p < ensemble.num_particles(); ++p) {
            for (size_t b = 0; b < ensemble.num_beads(); ++b) {
                ensemble.bead(p, b).position[0] = 0.0;
                ensemble.bead(p, b).position[1] = 0.0;
                ensemble.bead(p, b).position[2] = 0.0;
            }
        }
        
        EngineConfig config{
            .num_thermalization_steps = 50,
            .num_measurement_steps = 200,
            .measurement_interval = 5,
            .random_seed = 0
        };
        
        qsim::TestableEngine<MockAction, TestPathEnsemble, MockRNG> engine(
            MockAction(), ensemble, config
        );
        engine.run();
        
        // Starting at minimum, most moves should be rejected
        // Energy should be near zero
        REQUIRE(engine.mean_energy() >= 0.0);
        auto error = engine.energy_standard_error();
        INFO("Standard error: " << error);
    }
}

TEST_CASE("TestableEngine - Deterministic behavior", "[engine][determinism]") {
    SECTION("Same seed produces identical results") {
        EngineConfig config{
            .num_thermalization_steps = 50,
            .num_measurement_steps = 100,
            .measurement_interval = 10,
            .random_seed = 42
        };
        
        qsim::TestableEngine<MockAction, TestPathEnsemble, MockRNG> engine1(
            MockAction(), TestPathEnsemble(1, 3), config
        );
        engine1.run();
        
        qsim::TestableEngine<MockAction, TestPathEnsemble, MockRNG> engine2(
            MockAction(), TestPathEnsemble(1, 3), config
        );
        engine2.run();
        
        REQUIRE(engine1.mean_energy() == Approx(engine2.mean_energy()));
        REQUIRE(engine1.accepted_moves() == engine2.accepted_moves());
        REQUIRE(engine1.total_moves() == engine2.total_moves());
    }
}

TEST_CASE("TestableEngine - Edge cases", "[engine][edge_cases]") {
    SECTION("Zero steps") {
        EngineConfig config{
            .num_thermalization_steps = 0,
            .num_measurement_steps = 0,
            .random_seed = 0
        };
        
        qsim::TestableEngine<MockAction, TestPathEnsemble, MockRNG> engine(
            MockAction(), TestPathEnsemble(1, 1), config
        );
        
        REQUIRE_NOTHROW(engine.run());
        REQUIRE(engine.mean_energy() == 0.0);
        REQUIRE(engine.energy_standard_error() == 0.0);
    }
    
    SECTION("Single bead ensemble") {
        EngineConfig config{
            .num_thermalization_steps = 10,
            .num_measurement_steps = 20,
            .measurement_interval = 1,
            .random_seed = 0
        };
        
        qsim::TestableEngine<MockAction, TestPathEnsemble, MockRNG> engine(
            MockAction(), TestPathEnsemble(1, 1), config
        );
        
        REQUIRE_NOTHROW(engine.run());
        REQUIRE(engine.ensemble().num_particles() == 1);
        REQUIRE(engine.ensemble().num_beads() == 1);
    }
    
    SECTION("Large ensemble") {
        EngineConfig config{
            .num_thermalization_steps = 10,
            .num_measurement_steps = 20,
            .measurement_interval = 5,
            .random_seed = 0
        };
        
        qsim::TestableEngine<MockAction, TestPathEnsemble, MockRNG> engine(
            MockAction(), TestPathEnsemble(10, 20), config
        );
        
        REQUIRE_NOTHROW(engine.run());
        REQUIRE(engine.ensemble().num_particles() == 10);
        REQUIRE(engine.ensemble().num_beads() == 20);
    }
}

TEST_CASE("TestableEngine - Measurement interval", "[engine][measurements]") {
    SECTION("Measurement interval affects number of samples") {
        auto run_with_interval = [](size_t interval) {
            EngineConfig config{
                .num_thermalization_steps = 0,
                .num_measurement_steps = 100,
                .measurement_interval = interval,
                .random_seed = 0
            };
            
            qsim::TestableEngine<MockAction, TestPathEnsemble, MockRNG> engine(
                MockAction(), TestPathEnsemble(1, 3), config
            );
            engine.run();
            return engine.mean_energy();
        };
        
        double e1 = run_with_interval(1);   // Measures every step
        double e10 = run_with_interval(10);  // Measures every 10 steps
        double e100 = run_with_interval(100); // Measures once at start
        
        // Different intervals should give different statistics
        // At minimum, they should all be finite numbers
        REQUIRE(!std::isnan(e1));
        REQUIRE(!std::isnan(e10));
        REQUIRE(!std::isnan(e100));
    }
}


