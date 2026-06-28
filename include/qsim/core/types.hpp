// include/qsim/core/types.hpp
#ifndef QSIM_CORE_TYPES_HPP
#define QSIM_CORE_TYPES_HPP

#include <cstdint>
#include <array>
#include <vector>
#include <string>
#include <concepts>
#include <type_traits>
#include <compare>
#include<format>
#include <ostream>
#include <cmath>

namespace qsim {


//define integer types
using i8  = std::int8_t;
using i16 = std::int16_t;
using i32 = std::int32_t;
using i64 = std::int64_t;

using u8  = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

using usize = std::size_t;
using isize = std::ptrdiff_t;

using f32 = float;
using f64 = double;
using f128 = long double;  


// generic strong type to prevent unit mixing
template<typename T, typename Tag>
class StrongType {
public:
    using value_type = T;

    constexpr StrongType() noexcept : value_{} {}
    explicit constexpr StrongType(T value) noexcept : value_(value) {}

    // implicit conversion to undelying type
    constexpr operator T() const noexcept { return value_; }
    
    // accessors
    constexpr T value() const noexcept { return value_; }
    constexpr T& value() noexcept { return value_; }

    // comparison
    constexpr auto operator<=>(const StrongType& other) const = default;

    // type preserving arithmetic
    constexpr StrongType operator+(const StrongType& other) const {
        return StrongType(value_ + other.value_);
    }
    constexpr StrongType operator-(const StrongType& other) const {
        return StrongType(value_ - other.value_);
    }
    constexpr StrongType operator*(const StrongType& other) const {
        return StrongType(value_ * other.value_);
    }
    constexpr StrongType operator/(const StrongType& other) const {
        return StrongType(value_ / other.value_);
    }

    // scalar arithmetic 
    template<typename Scalar>
    constexpr StrongType operator*(Scalar scalar) const {
        return StrongType(value_ * scalar);
    }

    template<typename Scalar>
    friend constexpr StrongType operator*(Scalar scalar, const StrongType& v) {
        return StrongType(v.value_ * scalar);
    }

private:
    T value_;
};


// tag types 
struct MassTag {};
struct LengthTag {};
struct TimeTag {};
struct EnergyTag {};
struct TemperatureTag {};
struct ActionTag {};
struct FrequencyTag {};

// physical quantity 
using Mass       = StrongType<f64, MassTag>;
using Length     = StrongType<f64, LengthTag>;
using Time       = StrongType<f64, TimeTag>;
using Energy     = StrongType<f64, EnergyTag>;
using Temperature = StrongType<f64, TemperatureTag>;
using Action     = StrongType<f64, ActionTag>;
using Frequency  = StrongType<f64, FrequencyTag>;


//3d vector implementation 
template<typename T>
struct Vec3 {
    using value_type = T;
    T x, y, z;

    // constructors
    constexpr Vec3() noexcept : x(0), y(0), z(0) {}
    constexpr Vec3(T x, T y, T z) noexcept : x(x), y(y), z(z) {}
    explicit constexpr Vec3(T val) noexcept : x(val), y(val), z(val) {}

    // accessors
    constexpr T& operator[](usize i) noexcept {
        return (&x)[i];
    }
    constexpr const T& operator[](usize i) const noexcept {
        return (&x)[i];
    }

    // vector arithmetic
    constexpr Vec3 operator+(const Vec3& other) const {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }
    constexpr Vec3 operator-(const Vec3& other) const {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }
    constexpr Vec3 operator*(T scalar) const {
        return Vec3(x * scalar, y * scalar, z * scalar);
    }
    constexpr Vec3 operator/(T scalar) const {
        return Vec3(x / scalar, y / scalar, z / scalar);
    }

    // in place operators
    constexpr Vec3& operator+=(const Vec3& other) {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }
    constexpr Vec3& operator-=(const Vec3& other) {
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }
    constexpr Vec3& operator*=(T scalar) {
        x *= scalar; y *= scalar; z *= scalar;
        return *this;
    }

    // dot product
    constexpr T dot(const Vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    // cross product
    constexpr Vec3 cross(const Vec3& other) const {
        return Vec3(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }

    // norms
    constexpr T norm_sq() const {
        return x*x + y*y + z*z;
    }
    constexpr T norm() const {
        return std::sqrt(static_cast<long double>(norm_sq()));
    }

    // comparison
    constexpr bool operator==(const Vec3& other) const = default;
};

// common vector types
using Vec3f = Vec3<f32>;
using Vec3d = Vec3<f64>;


// strong type for particle indices
struct ParticleIndex {
    using value_type = i32;
    i32 id;

    constexpr ParticleIndex() noexcept : id(-1) {}
    explicit constexpr ParticleIndex(i32 id) noexcept : id(id) {}

    // implicit conversion to integer 
    constexpr operator i32() const noexcept { return id; }
    constexpr operator usize() const noexcept { return static_cast<usize>(id); }

    // explicit comparison 
    constexpr bool operator==(const ParticleIndex& other) const {
        return id == other.id;
    }
    constexpr bool operator!=(const ParticleIndex& other) const {
        return id != other.id;
    }
    constexpr bool operator<(const ParticleIndex& other) const {
        return id < other.id;
    }
    constexpr bool operator<=(const ParticleIndex& other) const {
        return id <= other.id;
    }
    constexpr bool operator>(const ParticleIndex& other) const {
        return id > other.id;
    }
    constexpr bool operator>=(const ParticleIndex& other) const {
        return id >= other.id;
    }

    // nullptr comparison 
    constexpr bool operator==(std::nullptr_t) const {
        return id == -1;
    }
    constexpr bool operator!=(std::nullptr_t) const {
        return id != -1;
    }





    // increment + decrement
    constexpr ParticleIndex& operator++() { ++id; return *this; }
    constexpr ParticleIndex operator++(int) { return ParticleIndex(id++); }
};

// strong type for bead time slices
struct BeadIndex {
    using value_type = i32;
    i32 index;

    constexpr BeadIndex() noexcept : index(-1) {}
    explicit constexpr BeadIndex(i32 idx) noexcept : index(idx) {}

    constexpr operator i32() const noexcept { return index; }
    constexpr operator usize() const noexcept { return static_cast<usize>(index); }

    constexpr auto operator<=>(const BeadIndex&) const = default;
};



// simulation parameters

struct SimulationParameters {
    f64 beta;           // inverse temperature (1/kt)
    i32 num_beads = 0;      // number of trotter slices
    i32 num_particles = 0;  // number of particles
    i32 num_dimensions = 3; // spatial dimensions (1, 2, or 3)
    f64 time_step;      // time step for propagation (beta / num_beads)
    
    // default  constructor
    SimulationParameters() = default;
    
    // constructor with validation
    constexpr SimulationParameters(
        f64 beta_in,
        i32 num_beads_in,
        i32 num_particles_in,
        i32 num_dimensions_in = 3
    ) : beta(beta_in),
        num_beads(num_beads_in),
        num_particles(num_particles_in),
        num_dimensions(num_dimensions_in),
        time_step(beta_in / static_cast<f64>(num_beads_in)) {
        
        // Validate parameters
        if (num_beads < 1) {
            throw std::invalid_argument("num_beads must be >= 1");
        }
        if (num_particles < 1) {
            throw std::invalid_argument("num_particles must be >= 1");
        }
        if (num_dimensions < 1 || num_dimensions > 3) {
            throw std::invalid_argument("num_dimensions must be 1, 2, or 3");
        }
        if (beta <= 0) {
            throw std::invalid_argument("beta must be positive");
        }
    }
    
    // comparison
    constexpr bool operator==(const SimulationParameters&) const = default;
};

// monte carlo statistics

struct MCStatistics {
    u64 total_sweeps = 0;
    u64 accepted_moves = 0;
    u64 rejected_moves = 0;
    f64 acceptance_rate = 0.0;
    f64 execution_time = 0.0;  // in seconds
    
    //reset
    void reset() {
        total_sweeps = 0;
        accepted_moves = 0;
        rejected_moves = 0;
        acceptance_rate = 0.0;
        execution_time = 0.0;
    }
    
    // record a move result
    void record_move(bool accepted) {
        ++total_sweeps;
        if (accepted) {
            ++accepted_moves;
        } else {
            ++rejected_moves;
        }
        update_acceptance_rate();
    }
    
    // update acceptance rate
    void update_acceptance_rate() {
        if (total_sweeps > 0) {
            acceptance_rate = static_cast<f64>(accepted_moves) / 
                             static_cast<f64>(total_sweeps);
        }
    }
};

// concepts for compile-time constraints

//  numeric types
template<typename T>
concept Numeric = std::is_arithmetic_v<T>;

//  floats
template<typename T>
concept Floating = std::is_floating_point_v<T>;

// real types
template<typename T>
concept Real = std::is_same_v<T, float> || 
               std::is_same_v<T, double> || 
               std::is_same_v<T, long double>;

// vector types
template<typename T>
concept Vector3 = requires(T v) {
    typename T::value_type;
    { v.x } -> std::convertible_to<typename T::value_type>;
    { v.y } -> std::convertible_to<typename T::value_type>;
    { v.z } -> std::convertible_to<typename T::value_type>;
    { v.dot(v) } -> std::convertible_to<typename T::value_type>;
    { v.norm() } -> std::convertible_to<typename T::value_type>;
};

// strong types with underlying value
template<typename T>
concept StrongTypeConcept = requires(T t) {
    typename T::value_type;
    { t.value() } -> std::convertible_to<typename T::value_type>;
};

} // end of namespace qsim

// custom formatters for strong types

//logging and debug 
namespace std {
    template<typename T, typename Tag>
    struct formatter<qsim::StrongType<T, Tag>> {
        constexpr auto parse(format_parse_context& ctx) {
            return ctx.begin();
        }
        
        auto format(const qsim::StrongType<T, Tag>& val, format_context& ctx) const {
            return format_to(ctx.out(), "{}", val.value());
        }
    };
}

#endif // QSIM_CORE_TYPES_HPP
