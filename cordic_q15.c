/*!
    \file
    \brief Optimized trigonometry functions.
    \authors Bliznets R.A. (r.bliznets@gmail.com)
    \version 1.0.0.0
    \date 11/11/2022
*/

#include <assert.h>
#include "mathlib_s3.h"
#include "complex_s3.h"
#include "sdkconfig.h"

// Lookup table for tangent values used in atan2 calculation (Q15 format)
static const q15 tan_array[16] = {4096, 2418, 1277, 648, 325, 163, 81, 41, 20, 10, 5, 2, 1, -16383, 1, 2};
// Internal implementation function for atan2 (potentially using specific hardware or optimized algorithm)
q15 atan2_q15_s3(q15 y, q15 x, const q15 *tan);
// Wrapper function for atan2 calculation using the lookup table
inline q15 atan2_q15(q15 y, q15 x)
{
    return atan2_q15_s3(y, x, tan_array);
}

// Internal implementation function for calculating arguments of a complex vector (PIE - Peripheral Input/Output Engine or specific hardware acceleration)
void arg_16_q15_pie(complex_q15 *in, q15 *out, uint32_t size, const q15 *tan);
// Wrapper function for calculating arguments of a complex vector with assertions
inline void arg_16_q15(complex_q15 *in, q15 *out, uint32_t size)
{
    // Ensure input vector is 16-byte aligned
    assert(((uint32_t)in % 16) == 0);
    // Ensure output vector is 16-byte aligned
    assert(((uint32_t)out % 16) == 0);
    // Ensure size is a multiple of 8
    assert((size % 8) == 0);
    // Ensure size is greater than 0
    assert(size > 0);

    arg_16_q15_pie(in, out, size, tan_array);
}

// Internal implementation function for sine and cosine calculation (potentially using specific hardware or optimized algorithm)
void sincos_q15_s3(q15 angle, q15 *cs, q15 *sn, const q15 *tan);
// Wrapper function for sine and cosine calculation using the lookup table
inline void sincos_q15(q15 angle, q15 *sn, q15 *cs)
{
    sincos_q15_s3(angle, cs, sn, tan_array);
}