// include/qsim/core/particle.hpp
#pragma once

#include <string>
#include <functional>
#include <cmath>
#include <compare>
#include <concepts>
#include <type_traits>

namespace qsim {

/**
 * @brief Represents a particle species with its quantum properties
 * 
 * This class defines the properties of a quantum particle species,
 * not individual particles in the simulation. It is immutable after
 * construction and serves as a type-safe way to specify particle
 * characteristics.
 */
class Particle {
public:
    /**
     * @brief Quantum statistics of the particle
     */
    enum class Statistics : uint8_t {
        Boson,      ///< Symmetric wavefunction, integer spin
        Fermion,    ///< Antisymmetric wavefunction, half-integer spin
        Boltzmann   ///< Classical (distinguishable) particles
    };

    /**
     * @brief Default constructor creates an invalid particle
     */
    Particle() = default;

    /**
     * @brief Construct a particle with specified properties
     * 
     * @param mass Mass in atomic units (electron mass = 1) AU
     * @param stats Quantum statistics
     * @param spin Spin quantum number (0 for bosons, 0.5 for fermions)
     * @param charge Charge in units of elementary charge
     */
    constexpr Particle(double mass, Statistics stats, 
                       double spin = 0.0, double charge = 0.0) noexcept
        : m_mass(mass)
        , m_spin(spin)
        , m_charge(charge)
        , m_statistics(stats)
    {}

    // Factory methods for common particles
    
    /**
     * @brief Create a bosonic particle
     * @param mass Mass in atomic units
     * @param charge Charge in units of e (default: 0)
     */
    static constexpr Particle boson(double mass, double charge = 0.0) noexcept {
        return Particle(mass, Statistics::Boson, 0.0, charge);
    }

    /**
     * @brief Create a fermionic particle
     * @param mass Mass in atomic units
     * @param spin Spin quantum number (default: 0.5)
     * @param charge Charge in units of e (default: 0)
     */
    static constexpr Particle fermion(double mass, double spin = 0.5, 
                                      double charge = 0.0) noexcept {
        return Particle(mass, Statistics::Fermion, spin, charge);
    }

    /**
     * @brief Create a classical (Boltzmann) particle
     * @param mass Mass in atomic units
     * @param charge Charge in units of e (default: 0)
     */
    static constexpr Particle boltzmann(double mass, double charge = 0.0) noexcept {
        return Particle(mass, Statistics::Boltzmann, 0.0, charge);
    }

    // Common species factories
    
    /**
     * @brief Helium-4 atom (boson)
     * @return Particle with He-4 properties
     */
    static constexpr Particle helium4() noexcept {
        return Particle(7294.299, Statistics::Boson, 0.0, 0.0);
    }

    /**
     * @brief Electron (fermion)
     * @return Particle with electron properties (mass in atomic units)
     */
    static constexpr Particle electron() noexcept {
        return Particle(1.0, Statistics::Fermion, 0.5, -1.0);
    }

    /**
     * @brief Proton (fermion)
     * @return Particle with proton properties (mass in atomic units)
     */
    static constexpr Particle proton() noexcept {
        return Particle(1836.152673, Statistics::Fermion, 0.5, 1.0);
    }

    /**
     * @brief Neutron (fermion)
     * @return Particle with neutron properties (mass in atomic units)
     */
    static constexpr Particle neutron() noexcept {
        return Particle(1838.683662, Statistics::Fermion, 0.5, 0.0);
    }

    // Accessors
    
    /**
     * @brief Get the mass of the particle
     * @return Mass in atomic units
     */
    constexpr double mass() const noexcept { return m_mass; }
    
    /**
     * @brief Get the spin quantum number
     * @return Spin (0 for bosons, 0.5 for fermions)
     */
    constexpr double spin() const noexcept { return m_spin; }
    
    /**
     * @brief Get the charge of the particle
     * @return Charge in units of elementary charge
     */
    constexpr double charge() const noexcept { return m_charge; }
    
    /**
     * @brief Get the quantum statistics
     * @return Statistics enum
     */
    constexpr Statistics statistics() const noexcept { return m_statistics; }

    /**
     * @brief Check if the particle is a boson
     * @return true if Statistics::Boson
     */
    constexpr bool is_boson() const noexcept { 
        return m_statistics == Statistics::Boson; 
    }

    /**
     * @brief Check if the particle is a fermion
     * @return true if Statistics::Fermion
     */
    constexpr bool is_fermion() const noexcept { 
        return m_statistics == Statistics::Fermion; 
    }

    /**
     * @brief Check if the particle is classical (Boltzmann)
     * @return true if Statistics::Boltzmann
     */
    constexpr bool is_boltzmann() const noexcept { 
        return m_statistics == Statistics::Boltzmann; 
    }

    /**
     * @brief Get string representation of statistics
     * @return "Boson", "Fermion", or "Boltzmann"
     */
    std::string statistics_string() const {
        switch (m_statistics) {
            case Statistics::Boson:    return "Boson";
            case Statistics::Fermion:  return "Fermion";
            case Statistics::Boltzmann: return "Boltzmann";
            default: return "Unknown";
        }
    }

    /**
     * @brief Check if the particle properties are physically valid
     * @return true if mass > 0 and statistics are consistent with spin
     */
    constexpr bool is_valid() const noexcept {
        if (m_mass <= 0.0) return false;
        
        // Bosons must have integer spin (0, 1, 2, ...)
        if (m_statistics == Statistics::Boson) {
            if (std::abs(m_spin - std::round(m_spin)) > 1e-12) return false;
            return true;
        }
        
        // Fermions must have half-integer spin (0.5, 1.5, ...)
        if (m_statistics == Statistics::Fermion) {
            double half_integer = m_spin - 0.5;
            if (std::abs(half_integer - std::round(half_integer)) > 1e-12) return false;
            return true;
        }
        
        // Boltzmann particles can have any spin
        return true;
    }

    // Comparison operators
    
    /**
     * @brief Equality comparison
     * @param other Particle to compare with
     * @return true if all properties are equal
     */
    constexpr bool operator==(const Particle& other) const noexcept {
        return m_mass == other.m_mass &&
               m_spin == other.m_spin &&
               m_charge == other.m_charge &&
               m_statistics == other.m_statistics;
    }

    /**
     * @brief Inequality comparison
     * @param other Particle to compare with
     * @return true if any property differs
     */
    constexpr bool operator!=(const Particle& other) const noexcept {
        return !(*this == other);
    }

    /**
     * @brief Three-way comparison for sorting
     * @param other Particle to compare with
     * @return std::strong_ordering based on mass, then spin, then statistics
     */
    constexpr std::partial_ordering operator<=>(const Particle& other) const noexcept {
        if (auto cmp = m_mass <=> other.m_mass; cmp != 0) return cmp;
        if (auto cmp = m_spin <=> other.m_spin; cmp != 0) return cmp;
        if (auto cmp = m_charge <=> other.m_charge; cmp != 0) return cmp;
        return m_statistics <=> other.m_statistics;
    }

private:
    double m_mass = 0.0;
    double m_spin = 0.0;
    double m_charge = 0.0;
    Statistics m_statistics = Statistics::Boltzmann;
};

/**
 * @brief Hash function for Particle class
 * 
 * Allows Particle to be used as a key in unordered containers
 */
struct ParticleHash {
    constexpr std::size_t operator()(const Particle& p) const noexcept {
        std::size_t h1 = std::hash<double>{}(p.mass());
        std::size_t h2 = std::hash<double>{}(p.spin());
        std::size_t h3 = std::hash<double>{}(p.charge());
        std::size_t h4 = std::hash<uint8_t>{}(
            static_cast<uint8_t>(p.statistics())
        );
        
        // Combine hashes
        return ((h1 ^ (h2 << 1)) >> 1) ^ 
               ((h3 ^ (h4 << 1)) >> 1);
    }
};

/**
 * @brief Concept for types that can be used as particles
 */
template<typename T>
concept ParticleConcept = requires(T p) {
    { p.mass() } -> std::convertible_to<double>;
    { p.spin() } -> std::convertible_to<double>;
    { p.charge() } -> std::convertible_to<double>;
    { p.statistics() } -> std::convertible_to<Particle::Statistics>;
    { p.is_valid() } -> std::convertible_to<bool>;
};

/**
 * @brief OStream operator for Particle
 */
inline std::ostream& operator<<(std::ostream& os, const Particle& p) {
    os << "Particle{mass=" << p.mass() 
       << ", spin=" << p.spin()
       << ", charge=" << p.charge()
       << ", stats=" << p.statistics_string() << "}";
    return os;
}

} // namespace qsim

