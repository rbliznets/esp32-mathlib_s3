#include "complex_s3.h"
#include "unity.h"

TEST_CASE("rotate_q15 periodic phases", "[mathlib]")
{
    constexpr uint32_t maxPeriod = 9;
    constexpr uint32_t samples = 12;
    const complex_q15 coefficients[maxPeriod] = {
        {16384, 0}, {0, 16384}, {-16384, 0}, {0, -16384}, {8192, 0},
        {0, 8192}, {-8192, 0}, {0, -8192}, {8192, 8192},
    };

    for (uint32_t period = 1; period <= maxPeriod; ++period)
    {
        // Zero rows past the logical table make a bad index observable without
        // letting a regressed implementation read outside the allocation.
        alignas(16) complex_q15 table[maxPeriod + 8][4] = {};
        for (uint32_t p = 0; p < period; ++p)
            for (uint32_t j = 0; j < 4; ++j)
                table[p][j] = coefficients[(p + j) % period];

        for (uint32_t phase = 0; phase < period; ++phase)
        {
            alignas(16) complex_q15 data[samples];
            complex_q15 expected[samples];
            for (uint32_t i = 0; i < samples; ++i)
            {
                data[i] = {static_cast<q15>(1024 + 128 * i), static_cast<q15>(-512 - 128 * static_cast<int>(i))};
                expected[i] = cmul_q15(data[i], coefficients[(phase + i) % period]);
            }
            rotate_q15(data, samples, table, period, phase);
            for (uint32_t i = 0; i < samples; ++i)
            {
                TEST_ASSERT_EQUAL_INT16(expected[i].re, data[i].re);
                TEST_ASSERT_EQUAL_INT16(expected[i].im, data[i].im);
            }
        }
    }
}
