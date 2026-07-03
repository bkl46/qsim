// include/qsim/rng/rng.hpp
#pragma once

#include <cstdint>
#include <type_traits>
#include <concepts>

namespace qsim {

/**
 * @brief Concept for Random Number Generator engines
 * 
 * Defines the minimal interface that any RNG engine must provide
 * to be used with the QSim framework.
 */
template<typename T>
concept RNGEngine = requires(T rng, std::array<uint64_t, 4> state) {
    // Must provide a default constructor
    typename T;
    
    // Must provide seed methods
    { rng.seed(uint64_t{}) } -> std::same_as<void>;
    
    // Must provide next() method for 64-bit random values
    { rng.next() } -> std::same_as<uint64_t>;
    
    // Must provide uniform distribution [0,1)
    { rng.uniform() } -> std::same_as<double>;
    
    // Must provide uniform distribution [0,1]
    { rng.uniform_inclusive() } -> std::same_as<double>;
    
    // Must provide Gaussian distribution (mean 0, std 1)
    { rng.gaussian() } -> std::same_as<double>;
    
    // Must provide state access
    { rng.get_state() } -> std::same_as<std::array<uint64_t, 4>>;
    
    // Must provide jump functionality for parallel streams
    { rng.jump() } -> std::same_as<void>;
    { rng.long_jump() } -> std::same_as<void>;
};

/**
 * @brief Base class for RNG engines
 * 
 * Provides common functionality and interface enforcement
 * for all RNG implementations in the QSim framework.
 */
class RNGBase {
public:
    virtual ~RNGBase() = default;
    
    /**
     * @brief Seed the RNG with a 64-bit value
     * @param seed The seed value
     */
    virtual void seed(uint64_t seed) = 0;
    
    /**
     * @brief Generate the next 64-bit random number
     * @return Uniformly distributed 64-bit unsigned integer
     */
    virtual uint64_t next() = 0;
    
    /**
     * @brief Generate a uniform random number in [0,1)
     * @return Double in range [0,1)
     */
    virtual double uniform() = 0;
    
    /**
     * @brief Generate a uniform random number in [0,1]
     * @return Double in range [0,1]
     */
    virtual double uniform_inclusive() = 0;
    
    /**
     * @brief Generate a Gaussian random number with mean 0, std 1
     * @return Standard normal random variable
     */
    virtual double gaussian() = 0;
    
    /**
     * @brief Get the current internal state
     * @return Array of 4 64-bit words representing the state
     */
    virtual std::array<uint64_t, 4> get_state() const = 0;
    
    /**
     * @brief Jump the RNG state forward by 2^128 steps
     */
    virtual void jump() = 0;
    
    /**
     * @brief Jump the RNG state forward by 2^192 steps
     */
    virtual void long_jump() = 0;
};

} // namespace qsim
