// include/qsim/path/path_ensemble.hpp

#pragma once

#include <vector>
#include <memory>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <fstream>
#include <cstring>
#include <stdexcept>
#include <concepts>
#include <type_traits>

#include "qsim/core/types.hpp"
#include "qsim/path/bead.hpp"
#include "qsim/path/path.hpp"
#include "qsim/core/constants.hpp"

namespace qsim {

class System;

/**
 * @brief Enumeration for memory layout strategies.
 */
enum class MemoryLayout {
    SoA,  // Structure of Arrays - better for SIMD
    AoS   // Array of Structures - easier to read
};

/**
 * @brief Concrete path ensemble class with optimized storage.
 */
class PathEnsemble {
public:
    // ========================================================================
    // Constructors
    // ========================================================================
    
    explicit PathEnsemble(const System& system, 
                          MemoryLayout layout = MemoryLayout::SoA)
        : system_ref_(&system)
        , layout_(layout) {
        reserve(system.get_num_particles(), system.get_num_beads(), system.get_dimensions());
    }
    
    PathEnsemble(size_t num_particles, uint32_t num_beads, 
                 uint32_t dimensions = 3,
                 MemoryLayout layout = MemoryLayout::SoA)
        : system_ref_(nullptr)
        , layout_(layout) {
        reserve(num_particles, num_beads, dimensions);
    }
    
    PathEnsemble() = default;
    
    ~PathEnsemble() = default;
    
    PathEnsemble(const PathEnsemble& other)
        : positions_(other.positions_)
        , topology_(other.topology_)
        , action_cache_(other.action_cache_)
        , num_particles_(other.num_particles_)
        , num_beads_(other.num_beads_)
        , dimensions_(other.dimensions_)
        , layout_(other.layout_)
        , bead_stride_(other.bead_stride_)
        , particle_stride_(other.particle_stride_)
        , dim_stride_(other.dim_stride_)
        , system_ref_(other.system_ref_) {}
    
    PathEnsemble(PathEnsemble&& other) noexcept
        : positions_(std::move(other.positions_))
        , topology_(std::move(other.topology_))
        , action_cache_(std::move(other.action_cache_))
        , num_particles_(other.num_particles_)
        , num_beads_(other.num_beads_)
        , dimensions_(other.dimensions_)
        , layout_(other.layout_)
        , bead_stride_(other.bead_stride_)
        , particle_stride_(other.particle_stride_)
        , dim_stride_(other.dim_stride_)
        , system_ref_(other.system_ref_) {
        other.num_particles_ = 0;
        other.num_beads_ = 0;
    }
    
    PathEnsemble& operator=(const PathEnsemble& other) {
        if (this != &other) {
            positions_ = other.positions_;
            topology_ = other.topology_;
            action_cache_ = other.action_cache_;
            num_particles_ = other.num_particles_;
            num_beads_ = other.num_beads_;
            dimensions_ = other.dimensions_;
            layout_ = other.layout_;
            bead_stride_ = other.bead_stride_;
            particle_stride_ = other.particle_stride_;
            dim_stride_ = other.dim_stride_;
            system_ref_ = other.system_ref_;
        }
        return *this;
    }
    
    PathEnsemble& operator=(PathEnsemble&& other) noexcept {
        if (this != &other) {
            positions_ = std::move(other.positions_);
            topology_ = std::move(other.topology_);
            action_cache_ = std::move(other.action_cache_);
            num_particles_ = other.num_particles_;
            num_beads_ = other.num_beads_;
            dimensions_ = other.dimensions_;
            layout_ = other.layout_;
            bead_stride_ = other.bead_stride_;
            particle_stride_ = other.particle_stride_;
            dim_stride_ = other.dim_stride_;
            system_ref_ = other.system_ref_;
            
            other.num_particles_ = 0;
            other.num_beads_ = 0;
        }
        return *this;
    }
    
    // ========================================================================
    // Memory Management
    // ========================================================================
    
    void reserve(size_t num_particles, uint32_t num_beads, uint32_t dimensions = 3) {
        num_particles_ = num_particles;
        num_beads_ = num_beads;
        dimensions_ = dimensions;
        
        const size_t total = num_particles * num_beads * dimensions;
        positions_.resize(total);
        action_cache_.resize(num_particles * num_beads);
        
        compute_strides();
    }
    
    void clear() {
        std::fill(positions_.begin(), positions_.end(), 0.0);
        std::fill(action_cache_.begin(), action_cache_.end(), 0.0);
    }
    
    void resize(size_t num_particles, uint32_t num_beads) {
        const size_t new_total = num_particles * num_beads * dimensions_;
        positions_.resize(new_total);
        action_cache_.resize(num_particles * num_beads);
        num_particles_ = num_particles;
        num_beads_ = num_beads;
        compute_strides();
    }
    
    // ========================================================================
    // Accessors
    // ========================================================================
    
    double get_bead(ParticleID particle_id, uint32_t bead_index, 
                    uint32_t dimension) const {
        return positions_[compute_index(particle_id, bead_index, dimension)];
    }
    
    double& get_bead(ParticleID particle_id, uint32_t bead_index, 
                     uint32_t dimension) {
        return positions_[compute_index(particle_id, bead_index, dimension)];
    }
    
    Vector3D get_bead_vector(ParticleID particle_id, uint32_t bead_index) const {
        Vector3D result{0.0, 0.0, 0.0};
        const size_t base = compute_index(particle_id, bead_index, 0);
        result.x = positions_[base];
        if (dimensions_ > 1) {
            result.y = positions_[base + (layout_ == MemoryLayout::SoA ? 1 : num_beads_)];
        }
        if (dimensions_ > 2) {
            result.z = positions_[base + (layout_ == MemoryLayout::SoA ? 2 : 2 * num_beads_)];
        }
        return result;
    }
    
    void set_bead(ParticleID particle_id, uint32_t bead_index, 
                  uint32_t dimension, double value) {
        positions_[compute_index(particle_id, bead_index, dimension)] = value;
    }
    
    void set_bead_vector(ParticleID particle_id, uint32_t bead_index, 
                         const Vector3D& position) {
        const size_t base = compute_index(particle_id, bead_index, 0);
        positions_[base] = position.x;
        if (dimensions_ > 1) {
            positions_[base + (layout_ == MemoryLayout::SoA ? 1 : num_beads_)] = position.y;
        }
        if (dimensions_ > 2) {
            positions_[base + (layout_ == MemoryLayout::SoA ? 2 : 2 * num_beads_)] = position.z;
        }
    }
    
    void set_positions(const std::vector<double>& positions) {
        if (positions.size() != positions_.size()) {
            throw std::invalid_argument("Position vector size mismatch");
        }
        std::copy(positions.begin(), positions.end(), positions_.begin());
    }
    
    // Raw data access
    double* data() { return positions_.data(); }
    const double* data() const { return positions_.data(); }
    size_t data_size() const { return positions_.size(); }
    
    // ========================================================================
    // Topology
    // ========================================================================
    
    const std::vector<ParticleID>& get_topology() const { return topology_; }
    
    void set_topology(const std::vector<ParticleID>& topology) {
        if (topology.size() != num_particles_) {
            throw std::invalid_argument("Topology size must match number of particles");
        }
        topology_ = topology;
    }
    
    void swap_paths(ParticleID particle_a, ParticleID particle_b) {
        if (particle_a >= topology_.size() || particle_b >= topology_.size()) {
            throw std::out_of_range("Invalid particle ID in swap_paths");
        }
        std::swap(topology_[particle_a], topology_[particle_b]);
    }
    
    // ========================================================================
    // Physics Computation
    // ========================================================================
    
    double compute_kinetic_action() const {
        if (!system_ref_) {
            throw std::runtime_error("System reference not set");
        }
        
        const double mass = system_ref_->get_particle(0).mass;
        const double tau = system_ref_->get_tau();
        const double prefactor = mass / (2.0 * constants::hbar * constants::hbar * tau);
        const double box_size = system_ref_->get_box_size();
        const bool periodic = (system_ref_->get_boundary() == BoundaryType::Periodic);
        
        double action = 0.0;
        
        for (size_t p = 0; p < num_particles_; ++p) {
            for (uint32_t b = 0; b < num_beads_; ++b) {
                const uint32_t next_b = (b + 1) % num_beads_;
                const Vector3D r1 = get_bead_vector(p, b);
                const Vector3D r2 = get_bead_vector(p, next_b);
                
                double dx = r2.x - r1.x;
                double dy = r2.y - r1.y;
                double dz = r2.z - r1.z;
                
                if (periodic) {
                    dx = dx - box_size * std::round(dx / box_size);
                    dy = dy - box_size * std::round(dy / box_size);
                    dz = dz - box_size * std::round(dz / box_size);
                }
                
                action += prefactor * (dx*dx + dy*dy + dz*dz);
            }
        }
        
        return action;
    }
    
    double compute_potential_action() const {
        if (!system_ref_) {
            throw std::runtime_error("System reference not set");
        }
        
        double potential = 0.0;
        const double tau = system_ref_->get_tau();
        
        for (uint32_t b = 0; b < num_beads_; ++b) {
            std::vector<double> bead_positions(num_particles_ * dimensions_);
            for (size_t p = 0; p < num_particles_; ++p) {
                for (uint32_t d = 0; d < dimensions_; ++d) {
                    bead_positions[p * dimensions_ + d] = get_bead(p, b, d);
                }
            }
            potential += system_ref_->evaluate_potential(bead_positions) * tau;
        }
        
        return potential;
    }
    
    double compute_total_action() const {
        return compute_kinetic_action() + compute_potential_action();
    }
    
    double compute_action_delta(ParticleID particle_id,
                                uint32_t bead_index,
                                const Vector3D& new_position) const {
        if (!system_ref_) {
            throw std::runtime_error("System reference not set");
        }
        
        const double mass = system_ref_->get_particle(particle_id).mass;
        const double tau = system_ref_->get_tau();
        const double prefactor = mass / (2.0 * constants::hbar * constants::hbar * tau);
        const double box_size = system_ref_->get_box_size();
        const bool periodic = (system_ref_->get_boundary() == BoundaryType::Periodic);
        
        const uint32_t prev_b = (bead_index == 0) ? num_beads_ - 1 : bead_index - 1;
        const uint32_t next_b = (bead_index + 1) % num_beads_;
        
        const Vector3D old_pos = get_bead_vector(particle_id, bead_index);
        const Vector3D prev_pos = get_bead_vector(particle_id, prev_b);
        const Vector3D next_pos = get_bead_vector(particle_id, next_b);
        
        auto kinetic_part = [&](const Vector3D& r1, const Vector3D& r2) -> double {
            double dx = r2.x - r1.x;
            double dy = r2.y - r1.y;
            double dz = r2.z - r1.z;
            if (periodic) {
                dx = dx - box_size * std::round(dx / box_size);
                dy = dy - box_size * std::round(dy / box_size);
                dz = dz - box_size * std::round(dz / box_size);
            }
            return prefactor * (dx*dx + dy*dy + dz*dz);
        };
        
        const double old_kinetic = kinetic_part(prev_pos, old_pos) + kinetic_part(old_pos, next_pos);
        const double new_kinetic = kinetic_part(prev_pos, new_position) + kinetic_part(new_position, next_pos);
        
        // Compute potential contribution
        std::vector<double> bead_positions(num_particles_ * dimensions_);
        for (size_t p = 0; p < num_particles_; ++p) {
            for (uint32_t d = 0; d < dimensions_; ++d) {
                if (p == particle_id) {
                    // Use old position initially
                    bead_positions[p * dimensions_ + d] = (d == 0) ? old_pos.x : 
                                                           (d == 1) ? old_pos.y : old_pos.z;
                } else {
                    bead_positions[p * dimensions_ + d] = get_bead(p, bead_index, d);
                }
            }
        }
        const double old_potential = system_ref_->evaluate_potential(bead_positions) * tau;
        
        // New potential
        for (uint32_t d = 0; d < dimensions_; ++d) {
            bead_positions[particle_id * dimensions_ + d] = (d == 0) ? new_position.x : 
                                                              (d == 1) ? new_position.y : new_position.z;
        }
        const double new_potential = system_ref_->evaluate_potential(bead_positions) * tau;
        
        return (new_kinetic - old_kinetic) + (new_potential - old_potential);
    }
    
    // ========================================================================
    // Modifiers
    // ========================================================================
    
    void update_bead(ParticleID particle_id, uint32_t bead_index,
                     const Vector3D& new_position, 
                     bool update_action_cache = true) {
        set_bead_vector(particle_id, bead_index, new_position);
        if (update_action_cache) {
            update_action_cache_impl(particle_id, bead_index);
        }
    }
    
    void update_beads(ParticleID particle_id,
                      const std::vector<uint32_t>& bead_indices,
                      const std::vector<Vector3D>& new_positions) {
        if (bead_indices.size() != new_positions.size()) {
            throw std::invalid_argument("Bead indices and positions size mismatch");
        }
        
        for (size_t i = 0; i < bead_indices.size(); ++i) {
            set_bead_vector(particle_id, bead_indices[i], new_positions[i]);
            update_action_cache_impl(particle_id, bead_indices[i]);
        }
    }
    
    void translate_all(const Vector3D& translation) {
        for (size_t p = 0; p < num_particles_; ++p) {
            for (uint32_t b = 0; b < num_beads_; ++b) {
                Vector3D pos = get_bead_vector(p, b);
                pos.x += translation.x;
                pos.y += translation.y;
                pos.z += translation.z;
                set_bead_vector(p, b, pos);
            }
        }
    }
    
    void apply_periodic_boundary(double box_size) {
        for (auto& pos : positions_) {
            pos = pos - box_size * std::round(pos / box_size);
        }
    }
    
    // ========================================================================
    // Estimators
    // ========================================================================
    
    Vector3D compute_path_center_of_mass(ParticleID particle_id) const {
        Vector3D com{0.0, 0.0, 0.0};
        for (uint32_t b = 0; b < num_beads_; ++b) {
            Vector3D pos = get_bead_vector(particle_id, b);
            com.x += pos.x;
            com.y += pos.y;
            com.z += pos.z;
        }
        com.x /= num_beads_;
        com.y /= num_beads_;
        com.z /= num_beads_;
        return com;
    }
    
    Vector3D compute_winding_number(ParticleID particle_id, double box_size) const {
        Vector3D winding{0.0, 0.0, 0.0};
        for (uint32_t b = 0; b < num_beads_; ++b) {
            const uint32_t next_b = (b + 1) % num_beads_;
            Vector3D r1 = get_bead_vector(particle_id, b);
            Vector3D r2 = get_bead_vector(particle_id, next_b);
            
            winding.x += std::round((r2.x - r1.x) / box_size);
            winding.y += std::round((r2.y - r1.y) / box_size);
            winding.z += std::round((r2.z - r1.z) / box_size);
        }
        return winding;
    }
    
    double compute_mean_squared_displacement(ParticleID particle_id, 
                                              uint32_t bead_separation) const {
        double msd = 0.0;
        const uint32_t half_beads = num_beads_ / 2;
        
        for (uint32_t b = 0; b < half_beads; ++b) {
            const uint32_t b2 = (b + bead_separation) % num_beads_;
            Vector3D r1 = get_bead_vector(particle_id, b);
            Vector3D r2 = get_bead_vector(particle_id, b2);
            
            double dx = r2.x - r1.x;
            double dy = r2.y - r1.y;
            double dz = r2.z - r1.z;
            
            msd += dx*dx + dy*dy + dz*dz;
        }
        
        return msd / half_beads;
    }
    
    // ========================================================================
    // Metadata
    // ========================================================================
    
    size_t get_num_particles() const { return num_particles_; }
    uint32_t get_num_beads() const { return num_beads_; }
    uint32_t get_dimensions() const { return dimensions_; }
    size_t get_total_beads() const { return num_particles_ * num_beads_; }
    bool empty() const { return positions_.empty(); }
    MemoryLayout get_memory_layout() const { return layout_; }
    void set_memory_layout(MemoryLayout layout) { layout_ = layout; }
    const System& get_system() const { return *system_ref_; }
    
    // ========================================================================
    // Serialization
    // ========================================================================
    
    void save_to_binary(const std::string& filename) const {
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Cannot open file for writing: " + filename);
        }
        
        uint64_t magic = 0x50494D43;  // "PIMC"
        file.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
        
        file.write(reinterpret_cast<const char*>(&num_particles_), sizeof(num_particles_));
        file.write(reinterpret_cast<const char*>(&num_beads_), sizeof(num_beads_));
        file.write(reinterpret_cast<const char*>(&dimensions_), sizeof(dimensions_));
        uint32_t layout_int = static_cast<uint32_t>(layout_);
        file.write(reinterpret_cast<const char*>(&layout_int), sizeof(layout_int));
        
        const size_t data_size = positions_.size();
        file.write(reinterpret_cast<const char*>(positions_.data()), data_size * sizeof(double));
        
        const size_t topo_size = topology_.size();
        file.write(reinterpret_cast<const char*>(topology_.data()), topo_size * sizeof(ParticleID));
    }
    
    static PathEnsemble load_from_binary(const std::string& filename, const System& system) {
        return load_from_binary_impl(filename, system);
    }
    
    static PathEnsemble load_from_binary_impl(const std::string& filename, const System& system) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Cannot open file for reading: " + filename);
        }
        
        uint64_t magic;
        file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        if (magic != 0x50494D43) {
            throw std::runtime_error("Invalid file format: " + filename);
        }
        
        size_t num_particles;
        uint32_t num_beads;
        uint32_t dimensions;
        uint32_t layout_int;
        
        file.read(reinterpret_cast<char*>(&num_particles), sizeof(num_particles));
        file.read(reinterpret_cast<char*>(&num_beads), sizeof(num_beads));
        file.read(reinterpret_cast<char*>(&dimensions), sizeof(dimensions));
        file.read(reinterpret_cast<char*>(&layout_int), sizeof(layout_int));
        
        MemoryLayout layout = static_cast<MemoryLayout>(layout_int);
        
        PathEnsemble ensemble(system, layout);
        ensemble.reserve(num_particles, num_beads, dimensions);
        
        const size_t data_size = num_particles * num_beads * dimensions;
        file.read(reinterpret_cast<char*>(ensemble.positions_.data()), 
                  data_size * sizeof(double));
        
        ensemble.topology_.resize(num_particles);
        file.read(reinterpret_cast<char*>(ensemble.topology_.data()), 
                  num_particles * sizeof(ParticleID));
        
        ensemble.num_particles_ = num_particles;
        ensemble.num_beads_ = num_beads;
        ensemble.dimensions_ = dimensions;
        ensemble.layout_ = layout;
        ensemble.compute_strides();
        
        return ensemble;
    }
    
    // ========================================================================
    // Comparison
    // ========================================================================
    
    bool operator==(const PathEnsemble& other) const {
        if (num_particles_ != other.num_particles_ ||
            num_beads_ != other.num_beads_ ||
            dimensions_ != other.dimensions_) {
            return false;
        }
        
        if (positions_ != other.positions_) return false;
        if (topology_ != other.topology_) return false;
        
        return true;
    }
    
    bool operator!=(const PathEnsemble& other) const { return !(*this == other); }
    
    bool is_valid() const {
        for (double pos : positions_) {
            if (!std::isfinite(pos)) return false;
        }
        
        if (topology_.size() != num_particles_) return false;
        if (positions_.size() != num_particles_ * num_beads_ * dimensions_) return false;
        
        return true;
    }

private:
    // ========================================================================
    // Private helper methods
    // ========================================================================
    
    size_t compute_index(ParticleID particle_id, uint32_t bead_index, 
                         uint32_t dimension) const {
        if (!is_valid_particle(particle_id) || !is_valid_bead(bead_index) || 
            dimension >= dimensions_) {
            throw std::out_of_range("Invalid index in PathEnsemble");
        }
        
        if (layout_ == MemoryLayout::SoA) {
            return (particle_id * num_beads_ + bead_index) * dimensions_ + dimension;
        } else {
            return (particle_id * dimensions_ + dimension) * num_beads_ + bead_index;
        }
    }
    
    void compute_strides() {
        if (layout_ == MemoryLayout::SoA) {
            dim_stride_ = 1;
            bead_stride_ = dimensions_;
            particle_stride_ = num_beads_ * dimensions_;
        } else {
            dim_stride_ = num_beads_;
            bead_stride_ = 1;
            particle_stride_ = dimensions_ * num_beads_;
        }
    }
    
    void update_action_cache_impl(ParticleID particle_id, uint32_t bead_index) {
        // Cache the action contribution for this bead
        // Implementation depends on what needs to be cached
        action_cache_[particle_id * num_beads_ + bead_index] = 0.0;
    }
    
    bool is_valid_bead(uint32_t bead_index) const {
        return bead_index < num_beads_;
    }
    
    bool is_valid_particle(ParticleID particle_id) const {
        return particle_id < num_particles_;
    }
    
    // ========================================================================
    // Member variables
    // ========================================================================
    
    std::vector<double> positions_;
    std::vector<ParticleID> topology_;
    std::vector<double> action_cache_;
    
    size_t num_particles_ = 0;
    uint32_t num_beads_ = 0;
    uint32_t dimensions_ = 3;
    MemoryLayout layout_ = MemoryLayout::SoA;
    
    size_t bead_stride_ = 0;
    size_t particle_stride_ = 0;
    size_t dim_stride_ = 0;
    
    const System* system_ref_ = nullptr;
};

// ========================================================================
// Free functions
// ========================================================================

inline void initialize_random(PathEnsemble& ensemble, double box_size) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-box_size/2.0, box_size/2.0);
    
    for (size_t i = 0; i < ensemble.data_size(); ++i) {
        ensemble.data()[i] = dist(gen);
    }
}

inline void initialize_straight_paths(PathEnsemble& ensemble) {
    // All beads at the same position (zero-temperature guess)
    // For each particle, get the centroid position and set all beads to it
    for (size_t p = 0; p < ensemble.get_num_particles(); ++p) {
        Vector3D pos = ensemble.get_bead_vector(p, 0);
        for (uint32_t b = 1; b < ensemble.get_num_beads(); ++b) {
            ensemble.set_bead_vector(p, b, pos);
        }
    }
}

inline double bead_distance(const PathEnsemble& ensemble, 
                            ParticleID particle_id,
                            uint32_t bead_i, uint32_t bead_j,
                            double box_size) {
    Vector3D r1 = ensemble.get_bead_vector(particle_id, bead_i);
    Vector3D r2 = ensemble.get_bead_vector(particle_id, bead_j);
    
    double dx = r2.x - r1.x;
    double dy = r2.y - r1.y;
    double dz = r2.z - r1.z;
    
    dx = dx - box_size * std::round(dx / box_size);
    dy = dy - box_size * std::round(dy / box_size);
    dz = dz - box_size * std::round(dz / box_size);
    
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

inline double kinetic_energy_between_beads(const Vector3D& r1, const Vector3D& r2,
                                           double mass, double tau) {
    double dx = r2.x - r1.x;
    double dy = r2.y - r1.y;
    double dz = r2.z - r1.z;
    
    return (mass / (2.0 * tau * tau)) * (dx*dx + dy*dy + dz*dz);
}

} // namespace qsim

