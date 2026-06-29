// include/qsim/path/path.hpp

#pragma once
#include<iostream>
#include <qsim/path/bead.hpp>
#include <qsim/core/types.hpp>
#include <vector>
#include <stdexcept>
#include <numeric>
#include <cmath>

namespace qsim {

/**
 * A ring polymer path for a single quantum particle
 * 
 * Stores P beads in a contiguous container,  path should satisfy cyclic boundary
 * conditions: the spring potential couples bead[i] to bead[(i+1)%P] for all i.
 * 
 * Key properties:
 *   - next(i) returns bead[(i+1) % P]
 *   - prev(i) returns bead[(i-1+P) % P]  
 *   - centroid() computes the imaginary-time average position
 * 
 */

template<typename T>
class Path {
public:
    using value_type = T;
    using bead_type = Bead<T>;
    using container_type = std::vector<bead_type>;
    using size_type = int;
    using iterator = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;


    //default constructor
    Path() noexcept = default;


    //construct path with P beads, all initialized to (0,0,0)
    explicit Path(int P) 
        : beads_(validate_size(P))
    {}


    //construct a path from a vector of of bead positions via copy
    explicit Path(const container_type& positions)
        : beads_(positions)
    {}


    //construct a path from a vector of of bead positions via move
    explicit Path(container_type&& positions) noexcept
        : beads_(std::move(positions))
    {}

    // default copy, move, and destructor
    Path(const Path&) = default;
    Path(Path&&) noexcept = default;
    Path& operator=(const Path&) = default;
    Path& operator=(Path&&) noexcept = default;
    ~Path() = default;

    // Capacity and access


    //num beads in path
    [[nodiscard]] constexpr size_type num_beads() const noexcept {
        return beads_.size();
    }


    //check empty
    [[nodiscard]] constexpr bool empty() const noexcept {
        return beads_.empty();
    }


    //access bead by time slice index (const)
    [[nodiscard]] const bead_type& bead(size_type i) const {
        return beads_.at(i);
    }


    //access bead by time slice index (mutable)
    [[nodiscard]] bead_type& bead(size_type i) {
        return beads_.at(i);
    }


    //set position for bead at time slice i
    void set_bead(size_type i, const bead_type& pos) {
        beads_.at(i) = pos;
    }

    // Cyclic navigation


    //get the next bead in imaginary time (wraps)
    //return beads[(i+1) % P] (const)
    [[nodiscard]] const bead_type& next(size_type i) const {
        if (i >= beads_.size()) {
            throw std::out_of_range("Path::next: index out of range");
        }
        return beads_[(i + 1) % beads_.size()];
    }



    //get the previous bead in imaginary time (wraps)
    //return beads[(i-1) % P] (const)
    [[nodiscard]] const bead_type& prev(size_type i) const {
        if (i >= beads_.size()) {
            throw std::out_of_range("Path::prev: index out of range");
        }
        return beads_[(i - 1 + beads_.size()) % beads_.size()];
    }

    // Physical quantities

    //compute the centroid, imagineary time average, of path (semi classical position)
    //  = (1/P) * sum_{i=0}^{P-1} bead[i]
    //  return vec representing centroid position
    [[nodiscard]] bead_type centroid() const {
        if (beads_.empty()) {
            throw std::runtime_error("Path::centroid: cannot compute centroid of empty path");
        }
        
        bead_type sum = std::accumulate(beads_.begin(), beads_.end(), bead_type{0});
        return sum / static_cast<T>(beads_.size());
    }

    // Iterator support for range-based for loops and algorithms

    [[nodiscard]] iterator begin() noexcept { return beads_.begin(); }
    [[nodiscard]] const_iterator begin() const noexcept { return beads_.begin(); }
    [[nodiscard]] const_iterator cbegin() const noexcept { return beads_.cbegin(); }
    
    [[nodiscard]] iterator end() noexcept { return beads_.end(); }
    [[nodiscard]] const_iterator end() const noexcept { return beads_.end(); }
    [[nodiscard]] const_iterator cend() const noexcept { return beads_.cend(); }

    // Direct container access for algorithms that need raw access


    //const reference to underlying bead container
    [[nodiscard]] const container_type& beads() const noexcept { return beads_; }

private:
    container_type beads_;  /// main container holding beads in time-slice order [0, P-1]
                            
                           
                          
    static int validate_size(int size) {
        if (size >= SIZE_MAX || size < 0) {
            throw std::invalid_argument("Invalid size parameter");
        }
        return size;
    }
};

// Type aliases for common floats

/// double-precision path 
using Path_d = Path<double>;

/// single-precision path 
using Path_f = Path<float>;

} // end of namespace qsim
