// include/qsim/path/bead.hpp

#pragma once

#include "../core/types.hpp"


namespace qsim {


/// bead in a path integral simulation, holds 3d position and time slice
/// 
/// Each bead represents the particle's position at a specific imaginary-time 
/// slice. The time-ordering is maintained by the Path container, thus bead 
/// itself carries no index information.

template<typename T>
using Bead = Vec3<T>;

/// Default bead type for double-precision simulations
using Bead_d = Bead<double>;

/// Bead type for single-precision simulations (GPU, memory-limited)
using Bead_f = Bead<float>;

} // namespace qsim
