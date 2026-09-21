/**
 * @file test_mathlib_s3_sum_pairs.cpp
 * @author Roman Bliznets (r.bliznets@gmail.com)
 * @brief Unity module test for sum_pairs_12_s16 (PIE sums of sample pairs on 12 adjacent shifts)
 * @version 0.0.0.1
 * @date 21.09.2026
 *
 * @copyright Copyright (c) 2026
 *
 */

#include <cstdint>
#include <cstring>
#include "unity.h"
#include "mathlib_s3.h"
#include "CTrace.h"
#include "esp_cpu.h"
#include "esp_random.h"

#define DATA_SIZE (1024)
#define MAX_PAIRS (64)

__attribute__((aligned(16))) static int16_t data[DATA_SIZE + 32];
static uint16_t pos[2 * MAX_PAIRS];

/// Reference: out[k] = sum_j (in[pos[2j+1] + k] - in[pos[2j] + k]).
static void sum_pairs_c(const int16_t *in, const uint16_t *p, uint32_t pairs, int32_t *out)
{
    for (int k = 0; k < 12; k++)
    {
        int32_t s = 0;
        for (uint32_t j = 0; j < pairs; j++)
            s += in[p[2 * j + 1] + k] - in[p[2 * j] + k];
        out[k] = s;
    }
}

TEST_CASE("sum_pairs_12_s16", "[mathlib]")
{
    __attribute__((aligned(16))) int32_t out[12];
    int32_t ref[12];

    for (int n = 0; n < 3000; n++)
    {
        // Full-scale extremes first, then random samples; random positions with every alignment.
        for (int i = 0; i < DATA_SIZE + 32; i++)
        {
            switch (n)
            {
            case 0:
                data[i] = INT16_MIN;
                break;
            case 1:
                data[i] = INT16_MAX;
                break;
            case 2:
                data[i] = (i & 1) ? INT16_MAX : INT16_MIN;
                break;
            default:
                data[i] = (int16_t)esp_random();
                break;
            }
        }
        uint32_t pairs = (n < 3) ? MAX_PAIRS : (esp_random() % (MAX_PAIRS + 1));
        for (uint32_t j = 0; j < 2 * pairs; j++)
        {
            if (n < 3)
                pos[j] = (j & 1) ? 1 : 0; // alternating extremes: differences of +-65535
            else
                pos[j] = (uint16_t)(esp_random() % (DATA_SIZE - 24 + 1));
        }
        std::memset(out, 0x55, sizeof(out));
        sum_pairs_12_s16(data, pos, pairs, out);
        sum_pairs_c(data, pos, pairs, ref);
        TEST_ASSERT_EQUAL_INT32_ARRAY(ref, out, 12);
    }

    // Timing: 31 pairs, as in the preamble correlation of PhaseDecoder.
    for (uint32_t j = 0; j < 62; j++)
        pos[j] = (uint16_t)(j * 10 + (j * 10) / 125);
    uint32_t c0 = esp_cpu_get_cycle_count();
    sum_pairs_12_s16(data, pos, 31, out);
    uint32_t c1 = esp_cpu_get_cycle_count();
    sum_pairs_c(data, pos, 31, ref);
    uint32_t c2 = esp_cpu_get_cycle_count();
    TEST_ASSERT_EQUAL_INT32_ARRAY(ref, out, 12);
    TRACE("sum_pairs_12_s16, 31 pairs, cycles", (int32_t)(c1 - c0), false);
    TRACE("C reference (12 shifts), cycles", (int32_t)(c2 - c1), false);
}
