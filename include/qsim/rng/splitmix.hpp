// include/qsim/rng/splitmix.hpp
#pragma once

#include <cstdint>
#include <array>

namespace qsim {

/**
 * @brief SplitMix64 Random Number Generator
 * 
 * Implementation of the SplitMix64 algorithm by Sebastiano Vigna.
 * This is a fast, high-quality 64-bit random number generator that is
 * primarily used for seeding other RNGs.
 * 
 * Reference: http://xorshift.di.unimi.it/splitmix64.c
 * 
 * Key features:
 * - 64-bit state
 * - Period of 2^64
 * - Excellent statistical properties
 * - Simple and fast
 * - Designed specifically for seeding generators like Xoshiro256
 */
class SplitMix64 {
public:
    /**
     * @brief Constructor with seed value
     * @param seed 64-bit seed value (default: 0)
     */
    explicit SplitMix64(uint64_t seed = 0) : state_(seed) {}
    
    /**
     * @brief Generate the next 64-bit random number
     * @return Uniformly distributed 64-bit unsigned integer
     * 
     * Uses the SplitMix64 algorithm with a mixed finalizer
     * to produce high-quality random numbers.
     */
    uint64_t next() {
        uint64_t z = (state_ += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    
    /**
     * @brief Generate a 64-bit random number in [0, 2^63-1]
     * @return Signed 64-bit integer in range [0, 2^63-1]
     */
    int64_t next_signed() {
        return static_cast<int64_t>(next() >> 1);
    }
    
    /**
     * @brief Generate a uniform random number in [0,1)
     * @return Double in range [0,1)
     */
    double uniform() {
        // Use 53 bits for double precision
        constexpr uint64_t MASK = (1ULL << 53) - 1;
        return static_cast<double>(next() & MASK) / (1ULL << 53);
    }
    
    /**
     * @brief Reset the generator with a new seed
     * @param seed New 64-bit seed value
     */
    void seed(uint64_t seed) {
        state_ = seed;
    }
    
    /**
     * @brief Get the current state
     * @return Current 64-bit state
     */
    uint64_t get_state() const {
        return state_;
    }
    
    /**
     * @brief Set the state directly
     * @param state New 64-bit state
     */
    void set_state(uint64_t state) {
        state_ = state;
    }
    
    /**
     * @brief Create a SplitMix64 instance seeded from the current state
     * @return New SplitMix64 instance
     * 
     * Useful for creating independent streams from a single seed.
     */
    SplitMix64 spawn() const {
        // Generate a new seed from the current state
        SplitMix64 temp(state_);
        uint64_t new_seed = temp.next();
        return SplitMix64(new_seed);
    }
    
private:
    uint64_t state_;
};

/**
 * @brief Utility class for seeding Xoshiro256 from SplitMix64
 * 
 * Provides static methods for generating seed states for
 * Xoshiro256 using the SplitMix64 algorithm.
 * 
 * Usage:
 *   auto state = Splitmix::seed(0x123456789ABCDEF0ULL);
 *   Xoshiro256 rng(state);
 */
class Splitmix {
public:
    /**
     * @brief Generate a seed state from a 64-bit seed
     * @param seed 64-bit seed value
     * @return Array of 4 64-bit words suitable for Xoshiro256
     * 
     * Uses SplitMix64 to expand a single 64-bit seed into
     * a 256-bit state (4 64-bit words).
     */
    static std::array<uint64_t, 4> seed(uint64_t seed) {
        std::array<uint64_t, 4> state;
        SplitMix64 sm(seed);
        
        for (size_t i = 0; i < 4; ++i) {
            state[i] = sm.next();
        }
        
        return state;
    }
    
    /**
     * @brief Generate a seed state from a seed sequence
     * @param seeds Array of seed values
     * @return Array of 4 64-bit words suitable for Xoshiro256
     * 
     * Uses SplitMix64 to expand an arbitrary number of seeds
     * into a 256-bit state.
     */
    template<typename Iterator>
    static std::array<uint64_t, 4> seed_from_sequence(Iterator begin, Iterator end) {
        std::array<uint64_t, 4> state = {0, 0, 0, 0};
        
        // Hash the sequence using SplitMix64
        for (auto it = begin; it != end; ++it) {
            SplitMix64 sm(*it);
            for (size_t i = 0; i < 4; ++i) {
                state[i] ^= sm.next();
            }
        }
        
        return state;
    }
    
    /**
     * @brief Generate a seed state from the system clock
     * @return Array of 4 64-bit words suitable for Xoshiro256
     * 
     * Uses time-based seeding for non-reproducible random sequences.
     * Useful for production simulations where reproducibility is
     * not required.
     */
    static std::array<uint64_t, 4> seed_from_time() {
        using namespace std::chrono;
        auto now = system_clock::now();
        uint64_t time_seed = duration_cast<nanoseconds>(now.time_since_epoch()).count();
        return seed(time_seed);
    }
    
    /**
     * @brief Generate a seed state from entropy sources
     * @return Array of 4 64-bit words suitable for Xoshiro256
     * 
     * Combines multiple entropy sources for robust seeding.
     * This is useful for production code where reproducibility
     * is not a concern.
     */
    static std::array<uint64_t, 4> seed_from_entropy() {
        // Combine time, address space randomization, and environment
        auto time_seed = duration_cast<nanoseconds>(
            std::chrono::system_clock::now().time_since_epoch()
        ).count();
        
        // Use stack address as a source of entropy
        uint64_t stack_seed = reinterpret_cast<uint64_t>(&stack_seed);
        
        // Use clock cycles for additional entropy
        uint64_t cycle_seed = static_cast<uint64_t>(__builtin_readcyclecounter());
        
        // Combine all sources
        uint64_t combined = time_seed ^ stack_seed ^ cycle_seed;
        
        return seed(combined);
    }
};

} // namespace qsim
