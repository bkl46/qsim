#pragma once

#include <concepts>
#include <type_traits>
#include "../path/path.hpp"
#include <span>
#include <cmath>
#include <memory>

namespace qsim {



/// @brief Concept defining the requirements for an action functional
/// @tparam A Action type
/// @tparam T Floating point type (typically double or float)
/// @tparam D Spatial dimension
template<typename A, typename T >
concept ActionConcept = requires(A a, const Path<T>& path, 
                                  //const PathEnsemble<T, D>& ensemble,
                                  size_t bead_idx, size_t particle_idx) {
    // Require floating point type for T
    requires std::is_floating_point_v<T>;
    
    // Core action evaluation methods
    { a.evaluate(path) } -> std::same_as<T>;
    { a.evaluate(ensemble) } -> std::same_as<T>;
    
    // Single bead action (for local updates)
    { a.evaluate_bead(path, bead_idx) } -> std::same_as<T>;
    { a.evaluate_bead(ensemble, particle_idx, bead_idx) } -> std::same_as<T>;
    
    // Action derivative for force-biased sampling
    { a.force(path, bead_idx) } -> std::same_as<std::array<T, D>>;
    { a.force(ensemble, particle_idx, bead_idx) } -> std::same_as<std::array<T, D>>;
    
    // Number of time slices (beads)
    { a.n_beads() } -> std::same_as<size_t>;
    
    // Clone support for polymorphic usage
    { a.clone() } -> std::same_as<std::unique_ptr<A>>;
    
    // Settable parameters
    { a.set_beta(1.0) } -> std::same_as<void>;
    { a.get_beta() } -> std::same_as<T>;
};

/// @brief Base class for all action functionals using CRTP
/// @tparam Derived The derived action type
/// @tparam T Floating point type for coordinates
/// @tparam D Spatial dimension (default 3)
template<typename Derived, typename T = double, size_t D = 3>
class ActionBase {
public:
    using value_type = T;
    using array_type = std::array<T, D>;
    static constexpr size_t dimension = D;

    ActionBase() : beta_(1.0), n_beads_(1) {}
    explicit ActionBase(T beta, size_t n_beads = 1) 
        : beta_(beta), n_beads_(n_beads) {}
    virtual ~ActionBase() = default;

    // Non-virtual interface (uses CRTP for dispatch)
    
    /// @brief Evaluate total action for a single path
    T evaluate(const Path<T, D>& path) const {
        return static_cast<const Derived*>(this)->evaluate_impl(path);
    }
    
    /// @brief Evaluate total action for path ensemble
    T evaluate(const PathEnsemble<T, D>& ensemble) const {
        return static_cast<const Derived*>(this)->evaluate_impl(ensemble);
    }
    
    /// @brief Evaluate action contribution from single bead
    T evaluate_bead(const Path<T, D>& path, size_t bead_idx) const {
        return static_cast<const Derived*>(this)->evaluate_bead_impl(path, bead_idx);
    }
    
    /// @brief Evaluate bead action in ensemble
    T evaluate_bead(const PathEnsemble<T, D>& ensemble, 
                    size_t particle_idx, size_t bead_idx) const {
        return static_cast<const Derived*>(this)->evaluate_bead_impl(
            ensemble, particle_idx, bead_idx);
    }
    
    /// @brief Compute force (negative gradient) at a bead
    array_type force(const Path<T, D>& path, size_t bead_idx) const {
        return static_cast<const Derived*>(this)->force_impl(path, bead_idx);
    }
    
    /// @brief Compute force for ensemble particle
    array_type force(const PathEnsemble<T, D>& ensemble,
                     size_t particle_idx, size_t bead_idx) const {
        return static_cast<const Derived*>(this)->force_impl(
            ensemble, particle_idx, bead_idx);
    }
    
    /// @brief Clone the action (virtual constructor)
    std::unique_ptr<Derived> clone() const {
        return std::make_unique<Derived>(
            static_cast<const Derived&>(*this));
    }

    // Accessors
    T get_beta() const noexcept { return beta_; }
    void set_beta(T beta) noexcept { beta_ = beta; }
    size_t n_beads() const noexcept { return n_beads_; }
    void set_n_beads(size_t n) noexcept { n_beads_ = n; }

    /// @brief Calculate time step tau = beta / n_beads
    T tau() const noexcept { return beta_ / static_cast<T>(n_beads_); }

protected:
    T beta_;         // Inverse temperature (1/kT)
    size_t n_beads_; // Number of time slices
};

/// @brief Tag types for different action contributions
namespace action_type {
    struct kinetic {};
    struct potential {};
    struct interaction {};
}

/// @brief Utility for combining action contributions
/// @tparam T Floating point type
template<typename T>
struct ActionResult {
    T total;
    T kinetic;
    T potential;
    T interaction;
    
    ActionResult& operator+=(const ActionResult& other) {
        total += other.total;
        kinetic += other.kinetic;
        potential += other.potential;
        interaction += other.interaction;
        return *this;
    }
};

} // namespace qsim
