#pragma once

#include "action.hpp"
#include "../path/path.hpp"
#include "../path/path_ensemble.hpp"
#include <array>
#include <cmath>
#include <algorithm>

namespace qsim {

/// @brief Primitive approximation for path integral action
/// 
/// The primitive approximation decomposes the density matrix as:
///   rho(R,R';tau) ≈ rho_free(R,R';tau) * exp(-tau * V(R))
/// 
/// For a path with M beads, the action is:
///   S = sum_{j=0}^{M-1} [ (m/2ħ²τ) * (R_j - R_{j+1})² + τ * V(R_j) ]
/// 
/// where indices are cyclic (R_M = R_0) for closed paths.
/// 
/// @tparam T Floating point type
/// @tparam D Spatial dimension
template<typename T = double, size_t D = 3>
class PrimitiveAction : public ActionBase<PrimitiveAction<T, D>, T, D> {
public:
    using Base = ActionBase<PrimitiveAction<T, D>, T, D>;
    using typename Base::value_type;
    using typename Base::array_type;
    using Base::dimension;
    using Base::beta_;
    using Base::n_beads_;

    /// @brief Construct primitive action with particle mass
    /// @param mass Particle mass (in units where ħ = 1)
    /// @param beta Inverse temperature
    /// @param n_beads Number of time slices
    explicit PrimitiveAction(T mass = 1.0, T beta = 1.0, size_t n_beads = 1)
        : Base(beta, n_beads), mass_(mass), 
          kinetic_prefactor_(mass / (2.0 * beta * beta)) {
        update_prefactors();
    }

    // --- Kinetic action ---
    
    /// @brief Compute total kinetic action for a single path
    T evaluate_impl(const Path<T, D>& path) const {
        T action = 0.0;
        const size_t M = path.n_beads();
        
        for (size_t j = 0; j < M; ++j) {
            const size_t jp1 = (j + 1) % M;  // Cyclic boundary
            action += distance_squared(path[j], path[jp1]);
        }
        
        return kinetic_prefactor_ * action;
    }
    
    /// @brief Compute total kinetic action for path ensemble
    T evaluate_impl(const PathEnsemble<T, D>& ensemble) const {
        T total_action = 0.0;
        for (size_t p = 0; p < ensemble.n_particles(); ++p) {
            total_action += evaluate_impl(ensemble[p]);
        }
        return total_action;
    }
    
    /// @brief Compute kinetic action contribution from a single bead
    T evaluate_bead_impl(const Path<T, D>& path, size_t bead_idx) const {
        const size_t M = path.n_beads();
        const size_t jm1 = (bead_idx + M - 1) % M;
        const size_t jp1 = (bead_idx + 1) % M;
        
        T action = distance_squared(path[jm1], path[bead_idx]) +
                   distance_squared(path[bead_idx], path[jp1]);
        
        return kinetic_prefactor_ * action;
    }
    
    /// @brief Compute kinetic action for a bead in ensemble
    T evaluate_bead_impl(const PathEnsemble<T, D>& ensemble,
                         size_t particle_idx, size_t bead_idx) const {
        return evaluate_bead_impl(ensemble[particle_idx], bead_idx);
    }
    
    // --- Forces (for force-biased Monte Carlo) ---
    
    /// @brief Compute kinetic force on a bead (negative gradient)
    array_type force_impl(const Path<T, D>& path, size_t bead_idx) const {
        const size_t M = path.n_beads();
        const size_t jm1 = (bead_idx + M - 1) % M;
        const size_t jp1 = (bead_idx + 1) % M;
        
        array_type force;
        const T coeff = 2.0 * kinetic_prefactor_;
        
        for (size_t d = 0; d < D; ++d) {
            force[d] = coeff * (path[jm1][d] - 2.0 * path[bead_idx][d] + path[jp1][d]);
        }
        
        return force;
    }
    
    /// @brief Compute kinetic force for ensemble particle
    array_type force_impl(const PathEnsemble<T, D>& ensemble,
                          size_t particle_idx, size_t bead_idx) const {
        return force_impl(ensemble[particle_idx], bead_idx);
    }
    
    // --- Configuration ---
    
    /// @brief Set particle mass and update prefactors
    void set_mass(T mass) noexcept {
        mass_ = mass;
        update_prefactors();
    }
    
    T get_mass() const noexcept { return mass_; }
    
    /// @brief Set inverse temperature and update prefactors
    void set_beta(T beta) noexcept {
        beta_ = beta;
        update_prefactors();
    }
    
    /// @brief Get the harmonic spring constant between beads
    T spring_constant() const noexcept { 
        return mass_ / (beta_ * tau()); 
    }
    
    /// @brief Standard deviation of free particle distribution between beads
    T free_particle_sigma() const noexcept {
        return std::sqrt(tau() / mass_);
    }

private:
    /// @brief Update cached prefactors when parameters change
    void update_prefactors() noexcept {
        const T tau = Base::tau();
        kinetic_prefactor_ = mass_ * static_cast<T>(n_beads_) / 
                            (2.0 * beta_);
    }
    
    /// @brief Compute squared Euclidean distance between two points
    static T distance_squared(const array_type& a, const array_type& b) noexcept {
        T dist2 = 0.0;
        for (size_t d = 0; d < D; ++d) {
            const T diff = a[d] - b[d];
            dist2 += diff * diff;
        }
        return dist2;
    }
    
    /// @brief Time step tau = beta / n_beads
    T tau() const noexcept { return beta_ / static_cast<T>(n_beads_); }

    T mass_;               // Particle mass
    T kinetic_prefactor_;  // Cached m*M/(2*beta) for kinetic term
};

// Verify concept satisfaction
static_assert(ActionConcept<PrimitiveAction<double, 3>, double, 3>);

} // namespace qsim
