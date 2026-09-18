/**
 * @file test_mathlib_s3_revbin.cpp
 * @author Roman Bliznets (r.bliznets@gmail.com)
 * @brief Unity module test for revbin_permute (PIE EE.BITREV for sizes 16..1024, C loop above)
 * @version 0.0.0.1
 * @date 18.09.2026
 *
 * @copyright Copyright (c) 2026
 *
 */

#include <cstdint>
#include <cstring>
#include <initializer_list>
#include "unity.h"
#include "fft_s3.h"
#include "CTrace.h"
#include "esp_random.h"
#include "esp_timer.h"

#define MAX_SIZE (4096)

__attribute__((aligned(16))) static complex_q15 data[MAX_SIZE];
__attribute__((aligned(16))) static complex_q15 ref[MAX_SIZE];

static uint32_t bitrev(uint32_t x, int bits)
{
    uint32_t r = 0;
    for (int i = 0; i < bits; i++)
        if (x & (1u << i))
            r |= 1u << (bits - 1 - i);
    return r;
}

/// Former C implementation of revbin_permute (kept for the timing comparison).
static void revbin_c(complex_q15 *d, uint32_t n)
{
    uint32_t *dt = (uint32_t *)d;
    uint32_t nh = n >> 1;
    uint32_t r = 0;
    uint32_t x = 1;
    uint32_t t;
    while (x < nh)
    {
        r = r + nh;
        t = dt[x];
        dt[x] = dt[r];
        dt[r] = t;
        x++;
        for (uint32_t m = n >> 1; (!((r ^= m) & m)); m >>= 1)
            ;
        if (r > x)
        {
            t = dt[x];
            dt[x] = dt[r];
            dt[r] = t;
            t = dt[n - 1 - x];
            dt[n - 1 - x] = dt[n - 1 - r];
            dt[n - 1 - r] = t;
        }
        x++;
    }
}

TEST_CASE("revbin_permute", "[mathlib_s3]")
{
    for (uint32_t n : {16u, 32u, 64u, 128u, 256u, 512u, 1024u, 2048u, 4096u})
    {
        int bits = fft_log2(n);
        for (uint32_t i = 0; i < n; i++)
        {
            uint32_t r = esp_random();
            std::memcpy(&data[i], &r, sizeof(r));
        }
        for (uint32_t i = 0; i < n; i++)
            ref[bitrev(i, bits)] = data[i];
        revbin_permute(data, n);
        TEST_ASSERT_EQUAL_MEMORY_MESSAGE(ref, data, n * sizeof(complex_q15), "revbin_permute");
        revbin_c(data, n); // the former C loop restores the natural order
        for (uint32_t i = 0; i < n; i++)
            TEST_ASSERT_EQUAL_UINT32(*(uint32_t *)&ref[bitrev(i, bits)], *(uint32_t *)&data[i]);
    }

    for (uint32_t n : {256u, 512u, 1024u})
    {
        int64_t t0 = esp_timer_get_time();
        for (int i = 0; i < 100; i++)
            revbin_permute(data, n);
        int64_t t1 = esp_timer_get_time();
        for (int i = 0; i < 100; i++)
            revbin_c(data, n);
        int64_t t2 = esp_timer_get_time();
        TRACE("revbin size", (int32_t)n, false);
        TRACE("  PIE, nsec", (int32_t)((t1 - t0) * 10), false);
        TRACE("  C, nsec", (int32_t)((t2 - t1) * 10), false);
    }
}
