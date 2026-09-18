/**
 * @file test_mathlib_s3_dot32.cpp
 * @author Roman Bliznets (r.bliznets@gmail.com)
 * @brief Unity module test for dot_product_16_16_32 (q15 dot product with a 32-bit result)
 * @version 0.0.0.1
 * @date 18.09.2026
 *
 * @copyright Copyright (c) 2026
 *
 */

#include <cstdint>
#include <cstring>
#include "unity.h"
#include "mathlib_s3.h"
#include "CTrace.h"
#include "esp_random.h"
#include "esp_timer.h"

#define MAX_SIZE (512)

__attribute__((aligned(16))) static q15 va[MAX_SIZE + 8];
__attribute__((aligned(16))) static q15 vb[MAX_SIZE + 8];

/// Reference: exact sum of products, arithmetic shift, saturation to int32.
static int32_t ref(const q15 *a, const q15 *b, uint32_t size, uint32_t shift)
{
    int64_t s = 0;
    for (uint32_t i = 0; i < size; i++)
        s += (int32_t)a[i] * (int32_t)b[i];
    s >>= shift;
    if (s > INT32_MAX)
        return INT32_MAX;
    if (s < INT32_MIN)
        return INT32_MIN;
    return (int32_t)s;
}

TEST_CASE("dot_product_16_16_32", "[mathlib_s3]")
{
    const uint32_t sizes[] = {8, 120, 512};
    const uint32_t shifts[] = {0, 1, 8, 15, 20};
    int32_t maxdiff = 0;
    for (int rep = 0; rep < 200; rep++)
    {
        for (uint32_t size : sizes)
        {
            // 512 products of 2^30 do not fit into the 40-bit accumulator
            if ((rep < 4) && (size > 256))
                continue;
            for (uint32_t i = 0; i < size + 8; i++)
            {
                uint32_t r = esp_random();
                if (rep < 4)
                {
                    // full scale: -32768 * -32768 and 32767 * -32768 products
                    va[i] = (rep & 1) ? INT16_MIN : INT16_MAX;
                    vb[i] = (rep & 2) ? INT16_MIN : INT16_MAX;
                }
                else
                {
                    va[i] = (q15)r;
                    vb[i] = (q15)(r >> 16);
                }
            }
            for (uint32_t shift : shifts)
            {
                int32_t d = dot_product_16_16_32(va, vb, size, shift) - ref(va, vb, size, shift);
                if (d < 0)
                    d = -d;
                if (d > maxdiff)
                    maxdiff = d;
            }
        }
    }
    TRACE("max |PIE - C|", maxdiff, false);

    for (uint32_t i = 0; i < 120; i++)
    {
        va[i] = (q15)esp_random();
        vb[i] = (q15)esp_random();
    }
    int64_t t0 = esp_timer_get_time();
    for (int i = 0; i < 1000; i++)
        dot_product_16_16_32(va, vb, 120, 8);
    TRACE("dot_product_16_16_32(120), nsec", (int32_t)(esp_timer_get_time() - t0), false);

    TEST_ASSERT_EQUAL_INT32(0, maxdiff);
}
