/*!
    \file
    \brief Optimized functions for complex numbers.
    \authors Bliznets R.A. (r.bliznets@gmail.com)
    \version 1.0.0.0
    \date 02/14/2023
*/

#include <assert.h>
#include "complex_s3.h"
#include "sdkconfig.h"

// Internal implementation function for magnitude calculation (PIE - Peripheral Input/Output Engine or specific hardware acceleration)
void magnitude_q15_pie(complex_q15 *in, q15 *out, uint32_t size);
// Wrapper function for magnitude calculation with assertions
inline void magnitude_q15(complex_q15 *in, q15 *out, uint32_t size)
{
    // Ensure input vector is 16-byte aligned
    assert(((uint32_t)in % 16) == 0);
    // Ensure output vector is 16-byte aligned
    assert(((uint32_t)out % 16) == 0);
    // Ensure size is a multiple of 8
    assert((size % 8) == 0);
    // Ensure size is greater than 0
    assert(size > 0);

    magnitude_q15_pie(in, out, size);
}

// Internal implementation function for complex multiplication (PIE)
complex_q15 cmul_q15_pie(complex_q15 x, complex_q15 y);
// Wrapper function for complex multiplication
inline complex_q15 cmul_q15(complex_q15 x, complex_q15 y)
{
    return cmul_q15_pie(x, y);
};

// Internal implementation function for multiplying a complex vector of size 10 by a scalar (PIE)
void cmul10_q15_pie(complex_q15 *in, complex_q15 *k, complex_q15 *out);
// Wrapper function for multiplying a complex vector of size 10 by a scalar with assertions
void cmul10_q15(complex_q15 *in, complex_q15 *k, complex_q15 *out)
{
    // Ensure input vector is 16-byte aligned
    assert(((uint32_t)in % 16) == 0);
    // Ensure scalar pointer is 4-byte aligned (for q15, which is int16_t, this ensures the complex_q15 struct is aligned)
    assert(((uint32_t)k % 4) == 0);
    // Ensure output vector is 16-byte aligned
    assert(((uint32_t)out % 16) == 0);

    cmul10_q15_pie(in, k, out);
}