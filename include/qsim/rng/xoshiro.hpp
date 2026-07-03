// include/qsim/rng/xoshiro.hpp
#pragma once

#include "rng.hpp"
#include <array>
#include <cstdint>
#include <cmath>
#include <limits>

namespace qsim {

/**
 * @brief Helper class for seeding Xoshiro256 from SplitMix64
 * 
 * This class provides static methods for generating seed states
 * for Xoshiro256 using the SplitMix64 algorithm.
 */
class Splitmix {
public:

    /**
     * @brief SplitMix64 seeding algorithm
     * @param seed 64-bit seed value
     * @return Array of 4 64-bit words representing the seeded state
     * 
     * Uses the SplitMix64 algorithm to expand a 64-bit seed
     * into a 256-bit state suitable for Xoshiro256.
     */
    static std::array<uint64_t, 4> splitmix64_seed(uint64_t seed) {
        std::array<uint64_t, 4> state;
        
        for (size_t i = 0; i < 4; ++i) {
            seed += 0x9e3779b97f4a7c15ULL;
            uint64_t z = seed;
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            z = z ^ (z >> 31);
            state[i] = z;
        }
        
        return state;
    }


    /**
     * @brief Generate a seed state from a 64-bit seed
     * @param seed 64-bit seed value
     * @return Array of 4 64-bit words suitable for Xoshiro256
     */
    static std::array<uint64_t, 4> seed(uint64_t seed) {
        return splitmix64_seed(seed);
    }
};





/**
 * @brief Xoshiro256** Random Number Generator
 * 
 * Implementation of the Xoshiro256** algorithm by David Blackman and Sebastiano Vigna.
 * This is a fast, high-quality RNG with good statistical properties and
 * jump functions for parallel computation.
 * 
 * Reference: https://prng.di.unimi.it/
 * 
 * Key features:
 * - 256-bit state (4 64-bit words)
 * - Period of 2^256 - 1
 * - Jump function for 2^128 streams
 * - Long jump for 2^192 streams
 * - Extremely fast (typically < 1-2 ns per generated number)
 * - Passes all TestU01 statistical tests
 */
class Xoshiro256 : public RNGBase {
public:
    /**
     * @brief Default constructor
     * 
     * Initializes with a fixed seed (0x123456789ABCDEF0ULL) for reproducibility.
     * Use seed() to initialize with a custom seed.
     */
    Xoshiro256() : Xoshiro256(0x123456789ABCDEF0ULL) {}
    
    /**
     * @brief Constructor with seed value
     * @param seed 64-bit seed value
     */
    explicit Xoshiro256(uint64_t seed) {
        this->seed(seed);
    }
    
    /**
     * @brief Constructor with explicit state
     * @param state Array of 4 64-bit words representing the state
     */
    explicit Xoshiro256(const std::array<uint64_t, 4>& state) 
        : s_(state) {}
    
    /**
     * @brief Copy constructor
     */
    Xoshiro256(const Xoshiro256& other) = default;
    
    /**
     * @brief Assignment operator
     */
    Xoshiro256& operator=(const Xoshiro256& other) = default;
    
    /**
     * @brief Seed the RNG with a 64-bit value
     * @param seed 64-bit seed value
     */
    void seed(uint64_t seed) override {
        // Use SplitMix64 to expand the 64-bit seed into 256-bit state
        s_ = Splitmix::splitmix64_seed(seed);
    }
    
    /**
     * @brief Generate the next 64-bit random number
     * @return Uniformly distributed 64-bit unsigned integer
     */
    uint64_t next() override {
        const uint64_t result = rotl(s_[1] * 5, 7) * 9;
        
        const uint64_t t = s_[1] << 17;
        
        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];
        
        s_[2] ^= t;
        s_[3] = rotl(s_[3], 45);
        
        return result;
    }
    
    /**
     * @brief Generate a uniform random number in [0,1)
     * @return Double in range [0,1)
     */
    double uniform() override {
        // Use the upper 53 bits for double precision
        constexpr uint64_t MASK = (1ULL << 53) - 1;
        return static_cast<double>(next() & MASK) / (1ULL << 53);
    }
    
    /**
     * @brief Generate a uniform random number in [0,1]
     * @return Double in range [0,1]
     */
    double uniform_inclusive() override {
        // Generate uniform in [0,1] with 53-bit precision
        constexpr uint64_t MASK = (1ULL << 53) - 1;
        uint64_t value = next() & MASK;
        // Add a small epsilon to ensure we can reach 1.0
        return static_cast<double>(value) / (1ULL << 53) + 
               (value == MASK ? 1.0 / (1ULL << 53) : 0.0);
    }
    
    /**
     * @brief Generate a Gaussian random number with mean 0, std 1
     * @return Standard normal random variable
     * 
     * Uses the Box-Muller transform with two uniforms per Gaussian.
     * This method is deterministic and thread-safe when called
     * from different RNG instances.
     */
    double gaussian() override {
        double u1, u2;
        do {
            u1 = uniform();  // (0,1]
        } while (u1 == 0.0);  // Avoid log(0)
        
        u2 = uniform();
        
        // Box-Muller transform
        const double r = std::sqrt(-2.0 * std::log(u1));
        const double theta = 2.0 * M_PI * u2;
        
        return r * std::cos(theta);
    }
    
    /**
     * @brief Get the current internal state
     * @return Array of 4 64-bit words representing the state
     */
    std::array<uint64_t, 4> get_state() const override {
        return s_;
    }
    
    /**
     * @brief Jump the RNG state forward by 2^128 steps
     * 
     * Used to create independent streams for parallel computations.
     * After jump(), the RNG produces a sequence that is independent
     * from the original sequence.
     * 
     * The jump polynomial is:
     * x^128 = 0x9e3779b185ebcaa7, 0xc2a1c1aa40c9b6c4, 
     *         0x8a28618cc0b6d19f, 0xb8b7de0de57c8d0c
     */
    void jump() override {
        constexpr std::array<uint64_t, 4> JUMP = {
            0x9e3779b185ebcaa7ULL,
            0xc2a1c1aa40c9b6c4ULL,
            0x8a28618cc0b6d19fULL,
            0xb8b7de0de57c8d0cULL
        };
        
        apply_jump(JUMP);
    }
    
    /**
     * @brief Jump the RNG state forward by 2^192 steps
     * 
     * Used to create independent streams for parallel computations.
     * Long jump provides even greater separation than jump().
     */
    void long_jump() override {
        constexpr std::array<uint64_t, 4> LONG_JUMP = {
            0xd978b1f9254ce4b4ULL,
            0x0783d7031e96980fULL,
            0x512aca38538f51b4ULL,
            0xc74e8fd4875231bdULL
        };
        
        apply_jump(LONG_JUMP);
    }
    
    /**
     * @brief Check if the RNG state is valid
     * @return true if the state is valid (not all zeros)
     */
    bool is_valid() const {
        return (s_[0] | s_[1] | s_[2] | s_[3]) != 0;
    }
    
private:
    std::array<uint64_t, 4> s_;  // Internal state
    
    /**
     * @brief Rotate-left operation
     * @param x 64-bit value to rotate
     * @param k Number of bits to rotate (mod 64)
     * @return Rotated value
     */
    static uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }
    
    /**
     * @brief Apply a jump polynomial to the state
     * @param jump_polynomial Array of 4 64-bit words representing the jump
     */
    void apply_jump(const std::array<uint64_t, 4>& jump_polynomial) {
        uint64_t s0 = 0;
        uint64_t s1 = 0;
        uint64_t s2 = 0;
        uint64_t s3 = 0;

        auto original = s_;
        
        for (size_t i = 0; i < 4; ++i) {
            for (int b = 0; b < 64; ++b) {
                if (jump_polynomial[i] & (1ULL << b)) {
                    s0 ^= original[0];
                    s1 ^= original[1];
                    s2 ^= original[2];
                    s3 ^= original[3];
                }
                // Advance the RNG by one step
                next();
            }
        }
        
        // Set the new state
        s_[0] = s0;
        s_[1] = s1;
        s_[2] = s2;
        s_[3] = s3;
    }
    


};


// Static assertion that Xoshiro256 satisfies the RNGEngine concept
static_assert(RNGEngine<Xoshiro256>);


} // namespace qsim
