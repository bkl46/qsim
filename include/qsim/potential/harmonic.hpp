// include/qsim/potential/harmonic.hpp

#pragma once
#include "../core/types.hpp"




namespace qsim::potential{



    class HarmonicOscillator{

        double m_omega; //frequency

        public:
        explicit HarmonicOscillator(double omega) : m_omega(omega) {}

        template <typename T> 
        double operator()(const Vec3<T>& x) const  {
            return 0.5 * m_omega * m_omega * (x.x*x.x + x.y*x.y + x.z*x.z);
        }


        template <typename  T>
        double gradient(const Vec3<T>& x) const {
            return m_omega * m_omega;
        }
    };



}
