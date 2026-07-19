// include/qsim/core/system.hpp

#pragma once

#include <vector>
#include <memory>
#include <string>
#include <unordered_map>
#include <random>
#include <fstream>
#include <numeric>
#include <type_traits>
#include <concepts>
#include <cstdint>
#include <algorithm>

#include "qsim/core/types.hpp"
#include "qsim/core/particle.hpp"
#include "qsim/core/constants.hpp"
#include "qsim/potential/potential.hpp"
#include "qsim/path/path_ensemble.hpp"

namespace qsim {

// Forward declaration
template<typename Derived>
class SystemBase;

/**
 * @brief CRTP base class for system definitions.
 */
template<typename Derived>
class SystemBase {
public:
    using ParticleContainer = std::vector<Particle>;
    using PotentialContainer = std::vector<std::unique_ptr<Potential>>;
    
    // Particle Management
    ParticleID add_particle(const Particle& particle) {
        return derived().add_particle_impl(particle);
    }
    
    std::vector<ParticleID> add_particles(const Particle& particle, size_t count) {
        return derived().add_particles_impl(particle, count);
    }
    
    const Particle& get_particle(ParticleID id) const {
        return derived().get_particle_impl(id);
    }
    
    const ParticleContainer& get_particles() const {
        return derived().get_particles_impl();
    }
    
    size_t get_num_particles() const {
        return derived().get_num_particles_impl();
    }
    
    bool are_identical(ParticleID a, ParticleID b) const {
        return derived().are_identical_impl(a, b);
    }
    
    // Potential Management
    void add_potential(std::unique_ptr<Potential> potential) {
        derived().add_potential_impl(std::move(potential));
    }
    
    void add_potential(std::unique_ptr<Potential> potential,
                       const std::vector<std::pair<ParticleID, ParticleID>>& pairs) {
        derived().add_potential_impl(std::move(potential), pairs);
    }
    
    const PotentialContainer& get_potentials() const {
        return derived().get_potentials_impl();
    }
    
    double evaluate_potential(const std::vector<double>& positions) const {
        return derived().evaluate_potential_impl(positions);
    }
    
    // Thermodynamic Parameters
    void set_temperature(double T) { derived().set_temperature_impl(T); }
    double get_temperature() const { return derived().get_temperature_impl(); }
    double get_beta() const { return derived().get_beta_impl(); }
    void set_num_beads(uint32_t N) { derived().set_num_beads_impl(N); }
    uint32_t get_num_beads() const { return derived().get_num_beads_impl(); }
    double get_tau() const { return derived().get_tau_impl(); }
    
    // Spatial Parameters
    void set_box_size(double L) { derived().set_box_size_impl(L); }
    double get_box_size() const { return derived().get_box_size_impl(); }
    void set_dimensions(uint32_t D) { derived().set_dimensions_impl(D); }
    uint32_t get_dimensions() const { return derived().get_dimensions_impl(); }
    void set_boundary(BoundaryType boundary) { derived().set_boundary_impl(boundary); }
    BoundaryType get_boundary() const { return derived().get_boundary_impl(); }
    
    // Factory Methods - Create PathEnsemble
    PathEnsemble create_initial_ensemble(
        InitializationType init_type = InitializationType::Random
    ) const {
        return derived().create_initial_ensemble_impl(init_type);
    }
    
    PathEnsemble create_ensemble_from_checkpoint(const std::string& filename) const {
        return derived().create_ensemble_from_checkpoint_impl(filename);
    }
    
    PathEnsemble create_ensemble_from_positions(
        const std::vector<double>& positions
    ) const {
        return derived().create_ensemble_from_positions_impl(positions);
    }
    
    PathEnsemble create_crystal_ensemble(double lattice_constant) const {
        return derived().create_crystal_ensemble_impl(lattice_constant);
    }
    
    // Serialization
    void save_to_json(const std::string& filename) const {
        derived().save_to_json_impl(filename);
    }
    
    static Derived load_from_json(const std::string& filename) {
        return Derived::load_from_json_impl(filename);
    }
    
    // Validation
    bool is_valid() const { return derived().is_valid_impl(); }
    std::string get_validation_errors() const { 
        return derived().get_validation_errors_impl(); 
    }

protected:
    Derived& derived() { return static_cast<Derived&>(*this); }
    const Derived& derived() const { return static_cast<const Derived&>(*this); }
    ~SystemBase() = default;
};

/**
 * @brief Concrete system class for PIMC simulations.
 */
class System : public SystemBase<System> {
public:
    System() = default;
    explicit System(double temperature) : temperature_(temperature) {
        beta_ = 1.0 / (constants::k_B * temperature_);
    }
    
    explicit System(double temperature, uint32_t num_beads) 
        : temperature_(temperature), num_beads_(num_beads) {
        beta_ = 1.0 / (constants::k_B * temperature_);
    }
    
    ~System() = default;
    
    System(const System&) = default;
    System(System&&) = default;
    System& operator=(const System&) = default;
    System& operator=(System&&) = default;
    
    // Implementation of SystemBase interface
    ParticleID add_particle_impl(const Particle& particle) {
        ParticleID id = next_particle_id_++;
        particles_.push_back(particle);
        particle_index_map_[id] = particles_.size() - 1;
        return id;
    }
    
    std::vector<ParticleID> add_particles_impl(const Particle& particle, size_t count) {
        std::vector<ParticleID> ids;
        ids.reserve(count);
        for (size_t i = 0; i < count; ++i) {
            ids.push_back(add_particle_impl(particle));
        }
        return ids;
    }
    
    const Particle& get_particle_impl(ParticleID id) const {
        auto it = particle_index_map_.find(id);
        if (it == particle_index_map_.end()) {
            throw std::out_of_range("Particle ID not found: " + std::to_string(id));
        }
        return particles_[it->second];
    }
    
    const ParticleContainer& get_particles_impl() const { return particles_; }
    size_t get_num_particles_impl() const { return particles_.size(); }
    
    bool are_identical_impl(ParticleID a, ParticleID b) const {
        const auto& pa = get_particle_impl(a);
        const auto& pb = get_particle_impl(b);
        return (pa.mass == pb.mass) && (pa.statistics == pb.statistics);
    }
    
    void add_potential_impl(std::unique_ptr<Potential> potential) {
        potentials_.push_back(std::move(potential));
        potential_pairs_.emplace_back();
    }
    
    void add_potential_impl(std::unique_ptr<Potential> potential,
                            const std::vector<std::pair<ParticleID, ParticleID>>& pairs) {
        potentials_.push_back(std::move(potential));
        potential_pairs_.push_back(pairs);
    }
    
    const PotentialContainer& get_potentials_impl() const { return potentials_; }
    
    double evaluate_potential_impl(const std::vector<double>& positions) const {
        double total_energy = 0.0;
        const size_t num_particles = particles_.size();
        const size_t dims = dimensions_;
        
        for (size_t p = 0; p < num_particles; ++p) {
            std::vector<double> pos(dims);
            for (size_t d = 0; d < dims; ++d) {
                pos[d] = positions[p * dims + d];
            }
            
            for (size_t pot_idx = 0; pot_idx < potentials_.size(); ++pot_idx) {
                const auto& pot = potentials_[pot_idx];
                const auto& pairs = potential_pairs_[pot_idx];
                
                if (pairs.empty()) {
                    total_energy += pot->evaluate(pos);
                } else {
                    for (const auto& pair : pairs) {
                        if (pair.first == p || pair.second == p) {
                            size_t other = (pair.first == p) ? pair.second : pair.first;
                            std::vector<double> other_pos(dims);
                            for (size_t d = 0; d < dims; ++d) {
                                other_pos[d] = positions[other * dims + d];
                            }
                            total_energy += pot->evaluate_pair(pos, other_pos);
                        }
                    }
                }
            }
        }
        
        return total_energy;
    }
    
    void set_temperature_impl(double T) {
        if (T <= 0.0) {
            throw std::invalid_argument("Temperature must be positive");
        }
        temperature_ = T;
        beta_ = 1.0 / (constants::k_B * T);
    }
    
    double get_temperature_impl() const { return temperature_; }
    double get_beta_impl() const { return beta_; }
    void set_num_beads_impl(uint32_t N) { num_beads_ = N; }
    uint32_t get_num_beads_impl() const { return num_beads_; }
    double get_tau_impl() const { return beta_ / static_cast<double>(num_beads_); }
    
    void set_box_size_impl(double L) { box_size_ = L; }
    double get_box_size_impl() const { return box_size_; }
    void set_dimensions_impl(uint32_t D) { dimensions_ = D; }
    uint32_t get_dimensions_impl() const { return dimensions_; }
    void set_boundary_impl(BoundaryType boundary) { boundary_ = boundary; }
    BoundaryType get_boundary_impl() const { return boundary_; }
    
    PathEnsemble create_initial_ensemble_impl(InitializationType init_type) const {
        PathEnsemble ensemble(*this);
        ensemble.reserve(particles_.size(), num_beads_, dimensions_);
        
        std::vector<double> positions;
        switch (init_type) {
            case InitializationType::Random:
                positions = generate_random_positions_impl();
                break;
            case InitializationType::Zero:
                positions.assign(particles_.size() * num_beads_ * dimensions_, 0.0);
                break;
            case InitializationType::Thermal:
            default:
                positions = generate_random_positions_impl();
                break;
        }
        
        ensemble.set_positions(positions);
        
        std::vector<ParticleID> topology(particles_.size());
        std::iota(topology.begin(), topology.end(), 0);
        ensemble.set_topology(topology);
        
        return ensemble;
    }
    
    PathEnsemble create_ensemble_from_positions_impl(
        const std::vector<double>& positions) const {
        
        if (!is_position_valid_impl(positions)) {
            throw std::invalid_argument("Invalid positions vector");
        }
        
        PathEnsemble ensemble(*this);
        ensemble.reserve(particles_.size(), num_beads_, dimensions_);
        ensemble.set_positions(positions);
        
        std::vector<ParticleID> topology(particles_.size());
        std::iota(topology.begin(), topology.end(), 0);
        ensemble.set_topology(topology);
        
        return ensemble;
    }
    
    PathEnsemble create_crystal_ensemble_impl(double lattice_constant) const {
        PathEnsemble ensemble(*this);
        ensemble.reserve(particles_.size(), num_beads_, dimensions_);
        
        std::vector<double> positions = generate_lattice_positions_impl(lattice_constant);
        ensemble.set_positions(positions);
        
        std::vector<ParticleID> topology(particles_.size());
        std::iota(topology.begin(), topology.end(), 0);
        ensemble.set_topology(topology);
        
        return ensemble;
    }
    
    PathEnsemble create_ensemble_from_checkpoint_impl(const std::string& filename) const {
        return PathEnsemble::load_from_binary_impl(filename, *this);
    }
    
    void save_to_json_impl(const std::string& filename) const {
        nlohmann::json j;
        
        j["particles"] = nlohmann::json::array();
        for (const auto& p : particles_) {
            j["particles"].push_back({
                {"name", p.name},
                {"mass", p.mass},
                {"statistics", static_cast<int>(p.statistics)},
                {"charge", p.charge}
            });
        }
        
        j["temperature"] = temperature_;
        j["num_beads"] = num_beads_;
        j["box_size"] = box_size_;
        j["dimensions"] = dimensions_;
        j["boundary"] = static_cast<int>(boundary_);
        
        std::ofstream file(filename);
        file << j.dump(4);
    }
    
    static System load_from_json_impl(const std::string& filename) {
        std::ifstream file(filename);
        nlohmann::json j;
        file >> j;
        
        System system;
        
        for (const auto& p : j["particles"]) {
            Statistics stats = static_cast<Statistics>(p["statistics"].get<int>());
            Particle particle(p["name"], p["mass"], stats, p["charge"]);
            system.add_particle_impl(particle);
        }
        
        system.set_temperature_impl(j["temperature"]);
        system.set_num_beads_impl(j["num_beads"]);
        system.set_box_size_impl(j["box_size"]);
        system.set_dimensions_impl(j["dimensions"]);
        system.set_boundary_impl(static_cast<BoundaryType>(j["boundary"].get<int>()));
        
        return system;
    }
    
    bool is_valid_impl() const {
        validation_errors_.clear();
        
        if (particles_.empty()) {
            validation_errors_ += "No particles defined. ";
        }
        
        if (temperature_ <= 0.0) {
            validation_errors_ += "Invalid temperature (must be > 0). ";
        }
        
        if (num_beads_ < 2) {
            validation_errors_ += "Number of beads must be at least 2. ";
        }
        
        if (box_size_ <= 0.0) {
            validation_errors_ += "Box size must be positive. ";
        }
        
        return validation_errors_.empty();
    }
    
    std::string get_validation_errors_impl() const {
        return validation_errors_;
    }

private:
    std::vector<double> generate_random_positions_impl() const {
        const size_t total_positions = particles_.size() * num_beads_ * dimensions_;
        std::vector<double> positions(total_positions);
        
        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(-box_size_/2.0, box_size_/2.0);
        
        for (auto& pos : positions) {
            pos = dist(gen);
        }
        
        return positions;
    }
    
    std::vector<double> generate_lattice_positions_impl(double lattice_constant) const {
        const size_t total_positions = particles_.size() * num_beads_ * dimensions_;
        std::vector<double> positions(total_positions);
        
        size_t particles_per_side = static_cast<size_t>(std::cbrt(particles_.size()));
        if (particles_per_side * particles_per_side * particles_per_side < particles_.size()) {
            particles_per_side++;
        }
        
        size_t idx = 0;
        for (size_t p = 0; p < particles_.size(); ++p) {
            size_t x = p % particles_per_side;
            size_t y = (p / particles_per_side) % particles_per_side;
            size_t z = p / (particles_per_side * particles_per_side);
            
            for (uint32_t b = 0; b < num_beads_; ++b) {
                positions[idx++] = (x - particles_per_side/2.0) * lattice_constant;
                positions[idx++] = (y - particles_per_side/2.0) * lattice_constant;
                if (dimensions_ > 2) {
                    positions[idx++] = (z - particles_per_side/2.0) * lattice_constant;
                }
            }
        }
        
        return positions;
    }
    
    bool is_position_valid_impl(const std::vector<double>& positions) const {
        const size_t expected = particles_.size() * num_beads_ * dimensions_;
        if (positions.size() != expected) {
            return false;
        }
        
        for (double pos : positions) {
            if (!std::isfinite(pos)) {
                return false;
            }
        }
        
        return true;
    }
    
    ParticleContainer particles_;
    std::unordered_map<ParticleID, size_t> particle_index_map_;
    ParticleID next_particle_id_ = 0;
    
    PotentialContainer potentials_;
    std::vector<std::vector<std::pair<ParticleID, ParticleID>>> potential_pairs_;
    
    double temperature_ = 300.0;
    double beta_ = 1.0 / (constants::k_B * 300.0);
    uint32_t num_beads_ = 64;
    
    double box_size_ = 10.0;
    uint32_t dimensions_ = 3;
    BoundaryType boundary_ = BoundaryType::Periodic;
    
    mutable std::string validation_errors_;
};

// Free functions
template<typename SystemType>
    requires std::derived_from<SystemType, SystemBase<SystemType>>
SystemType create_single_particle_system(const Particle& particle, double temperature) {
    SystemType system(temperature);
    system.add_particle(particle);
    return system;
}

template<typename SystemType>
    requires std::derived_from<SystemType, SystemBase<SystemType>>
SystemType create_fluid_system(const Particle& particle,
                               size_t num_particles,
                               double density,
                               double temperature) {
    SystemType system(temperature);
    system.add_particles(particle, num_particles);
    
    double volume = num_particles / density;
    double box_size = std::pow(volume, 1.0/3.0);
    system.set_box_size(box_size);
    
    return system;
}

} // namespace qsim

