
// tests/unit/test_rng_xoshiro.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
//#include <qsim/rng/xoshiro.hpp>
#include <qsim/rng/XoshiroCpp.hpp>
//#include <qsim/rng/splitmix.hpp>
#include <random>
#include <chrono>
#include <cmath>
#include <array>
#include <numeric>
#include <algorithm>
#include <iostream>

using namespace XoshiroCpp;
using Catch::Approx;

// Helper for statistical tolerance
constexpr double THREE_SIGMA = 3.0;
constexpr double FIVE_SIGMA = 5.0;
constexpr double TEN_SIGMA = 10.0;

TEST_CASE("Xoshiro256 RNG", "[rng][xoshiro]") {
    SECTION("Reproducibility - same seed yields same sequence") {
        const uint64_t seed = 0x123456789ABCDEF0ULL;
        Xoshiro256StarStar rng1(seed);
        Xoshiro256StarStar rng2(seed);
        
        std::array<uint64_t, 10> seq1, seq2;
        for (int i = 0; i < 10; ++i) {
            seq1[i] = rng1();
            seq2[i] = rng2();
        }
        
        REQUIRE(seq1 == seq2);
    }
    
    SECTION("Different seeds produce different sequences") {
        Xoshiro256StarStar rng1(0x123456789ABCDEF0ULL);
        Xoshiro256StarStar rng2(0xFEDCBA9876543210ULL);
        
        std::array<uint64_t, 10> seq1, seq2;
        for (int i = 0; i < 10; ++i) {
            seq1[i] = rng1();
            seq2[i] = rng2();
        }
        
        REQUIRE(seq1 != seq2);
    }
    
    //SECTION("Reseeding resets sequence") {
        //const uint64_t seed = 0x123456789ABCDEF0ULL;
        //Xoshiro256StarStar rng(seed);
        //
        //std::array<uint64_t, 10> seq1;
        //for (int i = 0; i < 10; ++i) {
            //seq1[i] = rng.next();
        //}
        //
        //// Reseed with same seed
        //rng.seed(seed);
        //std::array<uint64_t, 10> seq2;
        //for (int i = 0; i < 10; ++i) {
            //seq2[i] = rng.next();
        //}
        //
        //REQUIRE(seq1 == seq2);
    //}



}




TEST_CASE("Xoshiro256 uniform distribution", "[rng][xoshiro][statistics]") {
    SECTION("Basic range properties") {
        constexpr size_t NUM_SAMPLES = 100000000;
        Xoshiro256StarStar rng(0x123456789ABCDEF0ULL);
        
        uint64_t min_val = UINT64_MAX;
        uint64_t max_val = 0;
        double sum = 0.0;
        double sum_sq = 0.0;

        uint64_t value;
        
        for (size_t i = 0; i < NUM_SAMPLES; ++i) {
            value = rng();



            min_val = std::min(min_val, value);
            max_val = std::max(max_val, value);

            
            double u = DoubleFromBits(rng()); //rng.uniform();
            sum += u;
            sum_sq += u * u;
        }
        
        // Should cover nearly full range
        //REQUIRE(min_val <= 1000);
        //REQUIRE(max_val >= UINT64_MAX - 1000);
        //full range unlikely perhaps remake test
        
        // Test uniform [0,1) statistics
        double mean = sum / NUM_SAMPLES;
        double variance = (sum_sq / NUM_SAMPLES) - mean * mean;
        double expected_var = 1.0 / 12.0;
        double mean_sigma = std::sqrt(1.0 / 12.0 / NUM_SAMPLES);
        double var_sigma = std::sqrt(2.0 / (NUM_SAMPLES - 1)) * expected_var;
        
        REQUIRE(mean == Approx(0.5).margin(THREE_SIGMA * mean_sigma));
        REQUIRE(variance == Approx(expected_var).margin(THREE_SIGMA * var_sigma));
    }



    
    SECTION("Chi-square test for uniformity") {
        constexpr size_t NUM_SAMPLES = 10000;
        constexpr size_t NUM_BINS = 50;
        Xoshiro256StarStar rng(0x123456789ABCDEF0ULL);
        
        std::vector<size_t> bins(NUM_BINS, 0);
        for (size_t i = 0; i < NUM_SAMPLES; ++i) {
            double u = DoubleFromBits(rng());
            size_t bin = static_cast<size_t>(u * NUM_BINS);
            if (bin < NUM_BINS) {
                bins[bin]++;
            }
        }
        
        double expected = static_cast<double>(NUM_SAMPLES) / NUM_BINS;
        double chi2 = 0.0;
        for (size_t count : bins) {
            double diff = static_cast<double>(count) - expected;
            chi2 += diff * diff / expected;
        }
        
        // For 95% confidence with 49 DOF, chi2 ≈ 66.3
        // For 99% confidence with 49 DOF, chi2 ≈ 74.9
        REQUIRE(chi2 < 100.0);
    }


    
    SECTION("No significant autocorrelation") {
        constexpr size_t NUM_SAMPLES = 10000;
        Xoshiro256StarStar rng(0x123456789ABCDEF0ULL);
        
        std::vector<double> samples;
        samples.reserve(NUM_SAMPLES);
        for (size_t i = 0; i < NUM_SAMPLES; ++i) {
            samples.push_back(DoubleFromBits(rng()));
        }
        
        double mean = std::accumulate(samples.begin(), samples.end(), 0.0) / NUM_SAMPLES;
        double numerator = 0.0;
        double denominator = 0.0;
        
        for (size_t i = 0; i < NUM_SAMPLES - 1; ++i) {
            numerator += (samples[i] - mean) * (samples[i+1] - mean);
            denominator += (samples[i] - mean) * (samples[i] - mean);
        }
        
        double autocorr = numerator / denominator;
        double sigma = 1.0 / std::sqrt(NUM_SAMPLES);
        REQUIRE(autocorr == Approx(0.0).margin(FIVE_SIGMA * sigma));
    }

}



//write a gaussian function in Xoshiro

TEST_CASE("Xoshiro256 Gaussian distribution", "[rng][xoshiro][gaussian]") {
    SECTION("Basic statistical properties") {
        constexpr size_t NUM_SAMPLES = 1000000;
        Xoshiro256StarStar rng(0x123456789ABCDEF0ULL);
        
        double sum = 0.0;
        double sum_sq = 0.0;
        std::vector<double> samples;
        samples.reserve(NUM_SAMPLES);
        
        for (size_t i = 0; i < NUM_SAMPLES; ++i) {
            double g = gaussian(rng(), rng()); //rng.gaussian();
            REQUIRE(std::isfinite(g));
            sum += g;
            sum_sq += g * g;
            samples.push_back(g);
        }
        
        double mean = sum / NUM_SAMPLES;
        double variance = (sum_sq / NUM_SAMPLES) - mean * mean;
        
        double mean_sigma = std::sqrt(1.0 / NUM_SAMPLES);
        double var_sigma = std::sqrt(2.0 / (NUM_SAMPLES - 1));
        
        REQUIRE(mean == Approx(0.0).margin(FIVE_SIGMA * mean_sigma));
        REQUIRE(variance == Approx(1.0).margin(FIVE_SIGMA * var_sigma));
        
        // Test higher moments
        double sum_3 = 0.0, sum_4 = 0.0;
        for (double x : samples) {
            double z = x - mean;
            sum_3 += z * z * z;
            sum_4 += z * z * z * z;
        }
        double third_moment = sum_3 / NUM_SAMPLES;
        double fourth_moment = sum_4 / NUM_SAMPLES;
        
        double skew_sigma = std::sqrt(6.0 / NUM_SAMPLES);
        double kurt_sigma = std::sqrt(24.0 / NUM_SAMPLES);
        
        REQUIRE(third_moment == Approx(0.0).margin(TEN_SIGMA * skew_sigma));
        REQUIRE(fourth_moment == Approx(3.0).margin(TEN_SIGMA * kurt_sigma));
    }



    
    SECTION("Tail probabilities") {
        constexpr size_t NUM_SAMPLES = 100000000;
        Xoshiro256StarStar rng(0x123456789ABCDEF0ULL);
        
        size_t count_gt_2 = 0;
        size_t count_gt_3 = 0;
        size_t count_gt_4 = 0;
        size_t count_gt_5 = 0;
        
        for (size_t i = 0; i < NUM_SAMPLES; ++i) {
            double g = gaussian(rng(), rng()); //rng.gaussian();
            if (g > 2.0) count_gt_2++;
            if (g > 3.0) count_gt_3++;
            if (g > 4.0) count_gt_4++;
            if (g > 5.0) count_gt_5++;
        }
        
        auto check_tail = [&](double observed, double expected, size_t count) {
            printf("observed %d \n", observed);
            double sigma = std::sqrt(expected * (1 - expected) / NUM_SAMPLES);
            REQUIRE(observed == Approx(expected).margin(FIVE_SIGMA * sigma));
        };
        
        check_tail(static_cast<double>(count_gt_2) / NUM_SAMPLES, 0.02275, count_gt_2);
        check_tail(static_cast<double>(count_gt_3) / NUM_SAMPLES, 0.00135, count_gt_3);
        check_tail(static_cast<double>(count_gt_4) / NUM_SAMPLES, 0.0000317, count_gt_4);
        check_tail(static_cast<double>(count_gt_5) / NUM_SAMPLES, 0.000000287, count_gt_5);
    }


    
    SECTION("Symmetry around zero") {
        constexpr size_t NUM_SAMPLES = 100000;
        Xoshiro256StarStar rng(0x123456789ABCDEF0ULL);
        
        size_t count_positive = 0;
        size_t count_negative = 0;
        
        for (size_t i = 0; i < NUM_SAMPLES; ++i) {
            double g = gaussian(rng(), rng()); //rng.gaussian();
            if (g > 0) count_positive++;
            else if (g < 0) count_negative++;
        }
        
        double ratio = static_cast<double>(count_positive) / (count_positive + count_negative);
        double sigma = std::sqrt(0.25 / NUM_SAMPLES);
        REQUIRE(ratio == Approx(0.5).margin(FIVE_SIGMA * sigma));
    }


}



TEST_CASE("Xoshiro256 jump functionality", "[rng][xoshiro][jump]") {
    SECTION("jump() advances state deterministically") {
        const uint64_t seed = 0x123456789ABCDEF0ULL;
        Xoshiro256StarStar rng(seed);
        Xoshiro256StarStar rng2(seed);


        //same seed should have same state
        REQUIRE(rng == rng2);
        
        rng.jump();
        
        //different state after jump
        REQUIRE(rng != rng2);

        rng2.jump();

        //same state after jump
        REQUIRE(rng == rng2);
        
        std::array<uint64_t, 10> after_jump;
        for (int i = 0; i < 10; ++i) {
            after_jump[i] = rng();
        }
        
        //different state after use
        REQUIRE(rng != rng2);
        

    }
    

    SECTION("longJump() advances state deterministically") {
        const uint64_t seed = 0x123456789ABCDEF0ULL;
        Xoshiro256StarStar rng(seed);
        Xoshiro256StarStar rng2(seed);


        //same seed should have same state
        REQUIRE(rng == rng2);
        
        rng.longJump();
        
        //different state after jump
        REQUIRE(rng != rng2);

        rng2.longJump();

        //same state after jump
        REQUIRE(rng == rng2);
        
        std::array<uint64_t, 10> after_jump;
        for (int i = 0; i < 10; ++i) {
            after_jump[i] = rng();
        }
        
        //different state after use
        REQUIRE(rng != rng2);
        

    }
    


/*
    SECTION("long_jump() advances state by larger amount") {
        const uint64_t seed = 0x123456789ABCDEF0ULL;
        Xoshiro256StarStar rng(seed);
        
        std::array<uint64_t, 10> before_jump;
        for (int i = 0; i < 10; ++i) {
            before_jump[i] = rng.next();
        }
        
        auto state_before = rng.get_state();
        rng.long_jump();
        auto state_after = rng.get_state();
        
        REQUIRE(state_before != state_after);
        
        std::array<uint64_t, 10> after_jump;
        for (int i = 0; i < 10; ++i) {
            after_jump[i] = rng.next();
        }
        
        REQUIRE(before_jump != after_jump);
        
        // Long jump should be reproducible
        Xoshiro256 rng2(seed);
        rng2.long_jump();
        REQUIRE(rng.get_state() == rng2.get_state());
    }
    */

    
    SECTION("jump() and long_jump() produce different states") {
        const uint64_t seed = 0x123456789ABCDEF0ULL;
        Xoshiro256StarStar rng_jump(seed);
        Xoshiro256StarStar rng_long_jump(seed);

        REQUIRE(rng_jump == rng_long_jump);
        

        rng_jump.jump();
        rng_long_jump.longJump();
        
        REQUIRE(rng_jump != rng_long_jump);
    }


/*
    */
}



TEST_CASE("Xoshiro256 seeding", "[rng][xoshiro][seeding]") {
    SECTION("Splitmix seeding is reproducible") {
        auto seed1 = SplitMix64(0x123456789ABCDEF0ULL)();
        auto seed2 = SplitMix64(0x123456789ABCDEF0ULL)();
        
        REQUIRE(seed1 == seed2);
    }



    
    SECTION("Different Splitmix seeds produce different states") {
        auto state1 = SplitMix64(0x123456789ABCDEF0ULL)();
        auto state2 = SplitMix64(0xFEDCBA9876543210ULL)();
        
        REQUIRE(state1 != state2);
    }
    



    
    SECTION("Direct state initialization works") {

//

        SplitMix64 split = SplitMix64(0xFEDCBA9876543210ULL);

        std::uint64_t seed0 = split();
		std::uint64_t seed1 = split();
		std::uint64_t seed2 = split();
		std::uint64_t seed3 = split();
		Xoshiro256StarStar rng(
		{
		    seed0,
			seed1,
			seed2,
			seed3
		});

        std::array<std::uint64_t, 4> state = {seed0, seed1, seed2, seed3};
        REQUIRE(rng.serialize() == state);


        
        // Should generate valid random numbers
        for (int i = 0; i < 10; ++i) {
            uint64_t val = rng();
            REQUIRE(val != 0);  // Very unlikely to be zero
        }
    }

}





//curent work
TEST_CASE("Xoshiro256 serialization", "[rng][xoshiro][serialization]") {
    SECTION("State can be saved and restored") {
        const uint64_t seed = 0x123456789ABCDEF0ULL;
        Xoshiro256StarStar rng_original(seed);
        
        // Generate some random numbers
        for (int i = 0; i < 5; ++i) {
            rng_original();
        }
        
        auto saved_state = rng_original.serialize();
        
        // Create new RNG with saved state
        Xoshiro256StarStar rng_restored;
        rng_restored.deserialize(saved_state);
        
        // Generate sequences from both
        std::array<uint64_t, 20> seq1, seq2;
        for (int i = 0; i < 20; ++i) {
            seq1[i] = rng_original();
            seq2[i] = rng_restored();
        }
        
        // Sequences should be identical
        REQUIRE(seq1 == seq2);
    }
    

}



TEST_CASE("Xoshiro256 performance", "[rng][xoshiro][performance]") {
    SECTION("Generation speed is sufficient") {
        constexpr size_t NUM_ITERATIONS = 100000000;
        Xoshiro256StarStar rng(0x123456789ABCDEF0ULL);
        
        auto start = std::chrono::high_resolution_clock::now();
        
        volatile uint64_t sum = 0;
        for (size_t i = 0; i < NUM_ITERATIONS; ++i) {
            sum += rng();
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        // should generate at least 100m per second (10 ns per random)
        double ns_per_rand = static_cast<double>(duration.count()) * 1e6 / NUM_ITERATIONS;
        std::cout << "Xoshiro256 performance: " << ns_per_rand << " ns per random number\n";
        std::cout << "  Duration: " << duration << "\n";
        std::cout << " num iterations: " << NUM_ITERATIONS << "\n";
        REQUIRE(ns_per_rand < 10.0);
    }
}




TEST_CASE("Xoshiro256 edge cases", "[rng][xoshiro][edge]") {
    SECTION("Handles large number of samples") {
        constexpr size_t NUM_SAMPLES = 1000000;
        Xoshiro256StarStar rng(0x123456789ABCDEF0ULL);
        
        double sum_uniform = 0.0;
        double sum_gaussian = 0.0;
        
        for (size_t i = 0; i < NUM_SAMPLES; ++i) {
            sum_uniform += DoubleFromBits(rng());
            sum_gaussian += gaussian(rng(), rng());
        }
        
        REQUIRE(sum_uniform > NUM_SAMPLES * 0.49);
        REQUIRE(sum_uniform < NUM_SAMPLES * 0.51);
        REQUIRE(sum_gaussian > -1000.0);
        REQUIRE(sum_gaussian < 1000.0);
    }
    


    SECTION("All generated values are finite") {
        Xoshiro256StarStar rng(0x123456789ABCDEF0ULL);
        
        for (int i = 0; i < 1000; ++i) {
            REQUIRE(std::isfinite(DoubleFromBits(rng())));
            REQUIRE(std::isfinite(gaussian(rng(), rng())));
        }
    }
    
    SECTION("Handles zero seed gracefully") {
        Xoshiro256StarStar rng(0);
        
        // Should still generate good random numbers
        uint64_t sum = 0;
        for (int i = 0; i < 1000; ++i) {
            sum ^= rng();  // XOR to avoid overflow concerns
        }
        REQUIRE(sum != 0);  // Extremely unlikely to be zero
    }
}


TEST_CASE("Xoshiro256 stream generation", "[rng][xoshiro][stream]") {
    SECTION("Can generate random uint64_t values") {
        Xoshiro256StarStar rng(0x123456789ABCDEF0ULL);
        
        for (int i = 0; i < 100; ++i) {
            uint64_t val = rng();
            //  ensure it's valid
            REQUIRE(val >= 0);
        }
    }
    
    SECTION("Can generate random double in [0,1)") {
        Xoshiro256StarStar rng(0x123456789ABCDEF0ULL);
        
        for (int i = 0; i < 1000; ++i) {
            double u = DoubleFromBits(rng());
            REQUIRE(u >= 0.0);
            REQUIRE(u < 1.0);
        }
    }



}


