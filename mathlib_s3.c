/*!
    \file
    \brief Optimized DSP functions.
    \authors Bliznets R.A. (r.bliznets@gmail.com)
    \version 1.0.0.0
    \date 11/11/2022
*/

#include <assert.h>
#include "mathlib_s3.h"
#include "sdkconfig.h"

// Internal implementation function for vector copy (PIE - Peripheral Input/Output Engine or specific hardware acceleration)
void copy_pie(q15 *in, q15 *out, uint32_t size);
// Wrapper function for vector copy with assertions
inline void copy(q15 *in, q15 *out, uint32_t size)
{
    // Ensure output vector is 16-byte aligned
    assert(((uint32_t)out % 16) == 0);
    // Ensure size is a multiple of 8
    assert((size % 8) == 0);
    // Ensure size is greater than 0
    assert(size > 0);

    copy_pie(in, out, size);
}

// Internal implementation function for vector copy (16-byte aligned input) (PIE)
void copy_16_pie(q15 *in, q15 *out, uint32_t size);
// Wrapper function for vector copy (16-byte aligned input) with assertions
inline void copy_16(q15 *in, q15 *out, uint32_t size)
{
    // Ensure input vector is 16-byte aligned
    assert(((uint32_t)in % 16) == 0);
    // Ensure output vector is 16-byte aligned
    assert(((uint32_t)out % 16) == 0);
    // Ensure size is a multiple of 8
    assert((size % 8) == 0);
    // Ensure size is greater than 0
    assert(size > 0);

    copy_16_pie(in, out, size);
}

// Internal implementation function for multiplying vector by scalar (PIE)
void scaleVector_q15_pie(q15 *in, q15 *k, q15 *out, uint32_t size);
// Wrapper function for multiplying vector by scalar with assertions
inline void scaleVector(q15 *in, q15 *k, q15 *out, uint32_t size)
{
    // Ensure scalar pointer is 2-byte aligned (for q15)
    assert(((uint32_t)k % 2) == 0);
    // Ensure output vector is 16-byte aligned
    assert(((uint32_t)out % 16) == 0);
    // Ensure size is a multiple of 8
    assert((size % 8) == 0);
    // Ensure size is at least 16
    assert(size >= 16);

    scaleVector_q15_pie(in, k, out, size);
}

// Internal implementation function for multiplying vector by scalar (16-byte aligned input) (PIE)
void scaleVector_q15_16_pie(q15 *in, q15 *k, q15 *out, uint32_t size);
// Wrapper function for multiplying vector by scalar (16-byte aligned input) with assertions
inline void scaleVector_16(q15 *in, q15 *k, q15 *out, uint32_t size)
{
    // Ensure input vector is 16-byte aligned
    assert(((uint32_t)in % 16) == 0);
    // Ensure scalar pointer is 2-byte aligned (for q15)
    assert(((uint32_t)k % 2) == 0);
    // Ensure output vector is 16-byte aligned
    assert(((uint32_t)out % 16) == 0);
    // Ensure size is a multiple of 8
    assert((size % 8) == 0);
    // Ensure size is greater than 0
    assert(size > 0);

    scaleVector_q15_16_pie(in, k, out, size);
}

// Internal implementation function for shifting 32-bit vector to 16-bit vector (PIE)
void shrinkVector_16_pie(uint32_t *in, uint8_t shift, q15 *out, uint32_t size);
// Wrapper function for shifting 32-bit vector to 16-bit vector with assertions
inline void shrinkVector_16(uint32_t *in, uint8_t shift, q15 *out, uint32_t size)
{
    // Ensure input vector is 16-byte aligned
    assert(((uint32_t)in % 16) == 0);
    // Ensure shift amount is valid (less than 32)
    assert(shift < 32);
    // Ensure output vector is 16-byte aligned
    assert(((uint32_t)out % 16) == 0);
    // Ensure size is a multiple of 8
    assert((size % 8) == 0);
    // Ensure size is greater than 0
    assert(size > 0);

    shrinkVector_16_pie(in, shift, out, size);
}

// Internal implementation function for dot product (16-byte aligned vectors) (PIE)
q15 dot_product_q15_16_16_pie(q15 *in1, q15 *in2, uint32_t size);
// Wrapper function for dot product (16-byte aligned vectors) with assertions
inline q15 dot_product_16_16(q15 *in1, q15 *in2, uint32_t size)
{
    // Ensure first input vector is 16-byte aligned
    assert(((uint32_t)in1 % 16) == 0);
    // Ensure second input vector is 16-byte aligned
    assert(((uint32_t)in2 % 16) == 0);
    // Ensure size is a multiple of 8
    assert((size % 8) == 0);
    // Ensure size is greater than 0
    assert(size > 0);

    return dot_product_q15_16_16_pie(in1, in2, size);
}

// Internal implementation function for dot product (first vector unaligned) (PIE)
q15 dot_product_q15_1_16_pie(q15 *in1, q15 *in2, uint32_t size);
// Wrapper function for dot product (first vector unaligned) with assertions
inline q15 dot_product_1_16(q15 *in1, q15 *in2, uint32_t size)
{
    // Ensure second input vector is 16-byte aligned
    assert(((uint32_t)in2 % 16) == 0);
    // Ensure size is a multiple of 8
    assert((size % 8) == 0);
    // Ensure size is greater than 0
    assert(size > 0);

    return dot_product_q15_1_16_pie(in1, in2, size);
}

// FIR filter implementation using the unaligned dot product function.
void fir_1_16(q15 *in, q15 *k, uint32_t ksize, q15 *out, uint32_t size)
{
    // Ensure output size is greater than 0
    assert(size > 0);

    for (uint32_t i = 0; i < size; i++)
        // Calculate the dot product of the input window with the coefficients
        out[i] = dot_product_1_16(&in[i], k, ksize);
}

// Internal implementation function for adding two vectors (PIE)
void addVectors_q15_pie(q15 *in1, q15 *in2, q15 *out, uint32_t size);
// Wrapper function for adding two vectors with saturation and assertions
inline void addVectors_q15(q15 *in1, q15 *in2, q15 *out, uint32_t size)
{
    // Ensure first input vector is 16-byte aligned
    assert(((uint32_t)in1 % 16) == 0);
    // Ensure second input vector is 16-byte aligned
    assert(((uint32_t)in2 % 16) == 0);
    // Ensure output vector is 16-byte aligned
    assert(((uint32_t)out % 16) == 0);
    // Ensure size is a multiple of 8
    assert((size % 8) == 0);
    // Ensure size is greater than 0
    assert(size > 0);

    addVectors_q15_pie(in1, in2, out, size);
}

// Internal implementation function for normalizing vector q15 (PIE)
int16_t normalize_q15_pie(q15 *in, q15 *out, uint32_t size);
// Wrapper function for normalizing vector q15 with assertions
inline int16_t normalize_q15(q15 *in, q15 *out, uint32_t size)
{
    // Ensure input vector is 16-byte aligned
    assert(((uint32_t)in % 16) == 0);
    // Ensure output vector is 16-byte aligned
    assert(((uint32_t)out % 16) == 0);
    // Ensure size is a multiple of 8
    assert((size % 8) == 0);
    // Ensure size is greater than 0
    assert(size > 0);

    return normalize_q15_pie(in, out, size);
}

// Internal implementation function for normalizing vector q14 (PIE)
int16_t normalize_q14_pie(q15 *in, q15 *out, uint32_t size);
// Wrapper function for normalizing vector q14 with assertions
inline int16_t normalize_q14(q15 *in, q15 *out, uint32_t size)
{
    // Ensure input vector is 16-byte aligned
    assert(((uint32_t)in % 16) == 0);
    // Ensure output vector is 16-byte aligned
    assert(((uint32_t)out % 16) == 0);
    // Ensure size is a multiple of 8
    assert((size % 8) == 0);
    // Ensure size is greater than 0
    assert(size > 0);

    return normalize_q14_pie(in, out, size);
}

// Internal implementation function for FIR filter (PIE)
void fir_16_16_q15_pie(q15 *in, q15 *k, uint32_t ksize, q15 *out, uint32_t size);
// Wrapper function for FIR filter with assertions
inline void fir_16_16_q15(q15 *in, q15 *k, uint32_t ksize, q15 *out, uint32_t size)
{
    // Ensure input data vector is 16-byte aligned
    assert(((uint32_t)in % 16) == 0);
    // Ensure coefficient vector is 16-byte aligned
    assert(((uint32_t)k % 16) == 0);
    // Ensure output vector is 16-byte aligned
    assert(((uint32_t)out % 16) == 0);
    // Ensure coefficient size is a multiple of 8
    assert((ksize % 8) == 0);
    // Ensure coefficient size is greater than 0
    assert(ksize > 0);
    // Ensure output size is a multiple of 8
    assert((size % 8) == 0);
    // Ensure output size is greater than 8
    assert(size > 8);

    fir_16_16_q15_pie(in, k, ksize, out, size);
}