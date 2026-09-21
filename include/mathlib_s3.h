/*!
    \file
    \brief Optimized DSP functions.
    \authors Bliznets R.A. (r.bliznets@gmail.com)
    \version 1.0.0.0
    \date 11/11/2022
*/

#pragma once
#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

#define M_PI 3.14159265358979323846

    typedef int16_t q15; ///< fixed point S.15
#define toQ15(f) ((q15)(f * INT16_MAX))
#define toFloat(q) ((((int16_t)q)) / ((float)INT16_MAX))

#define toQ15Angle(f) ((q15)((f / (2 * M_PI)) * INT16_MAX))
#define toFloatAngle(q) ((q / float(INT16_MAX)) * 2 * M_PI)

#define to16SIZE(size) ((((size) + 15) / 16) * 16)
#define to8SIZE(size) ((((size) + 7) / 8) * 8)
#define to4SIZE(size) ((((size) + 3) / 4) * 4)

    /// Copy vector.
    /*!
        out=in
        \param[in] in vector.
        \param[out] out output vector (16 bytes aligned).
        \param[in] size vector size (multiple of 8, >= 8).
    */
    void copy(q15 *in, q15 *out, uint32_t size);
    /// Copy vector.
    /*!
        out=in
        \param[in] in vector (16 bytes aligned).
        \param[out] out output vector (16 bytes aligned).
        \param[in] size vector size (multiple of 8, >= 8).
    */
    void copy_16(q15 *in, q15 *out, uint32_t size);

    /// Multiply vector by scalar.
    /*!
        out=in * k
        \param[in] in vector (16 bytes aligned).
        \param[in] k pointer to scalar (2 bytes aligned).
        \param[out] out output vector (16 bytes aligned).
        \param[in] size vector size (multiple of 8, >= 8).
    */
    void scaleVector_16(q15 *in, q15 *k, q15 *out, uint32_t size);

    /// Shift 32bit vector to 16bit vector.
    /*!
        out=in >> shift
        \param[in] in 32bit vector (16 bytes aligned).
        \param[in] shift right shift (1..31).
        \param[out] out output 15bit vector (16 bytes aligned).
        \param[in] size vector size (multiple of 8, >= 8).
    */
    void shrinkVector_16(uint32_t *in, uint8_t shift, q15 *out, uint32_t size);

    /// Maximum absolute value of 32bit vector shifted to 16bit.
    /*!
        Bits 30..15 of each word are taken, bit 31 is ignored (AK2401 IQ words: the 24bit sample
        occupies bits 30..7, bit 31 is the tail of the previous I2S slot).
        \param[in] in 32bit vector (16 bytes aligned).
        \param[in] size vector size in words (multiple of 8, >= 8).
        \return max abs of (int16_t)(in[i] >> 15).
    */
    uint32_t maxAbsVector_16(uint32_t *in, uint32_t size);

    /// Multiply vector by scalar.
    /*!
        out=in * k
        \param[in] in vector.
        \param[in] k pointer to scalar (2 bytes aligned).
        \param[out] out output vector (16 bytes aligned).
        \param[in] size vector size (multiple of 8, >= 16).
    */
    void scaleVector(q15 *in, q15 *k, q15 *out, uint32_t size);

    /// Dot product of vectors q15.
    /*!
        \param[in] in1 vector (16 bytes aligned).
        \param[in] in2 vector (16 bytes aligned).
        \param[in] size vector size (multiple of 8, >= 8).
        \return result
    */
    q15 dot_product_16_16(q15 *in1, q15 *in2, uint32_t size);
    /// Dot product of vectors q15 with a 32-bit result.
    /*!
        The exact sum of products is accumulated in the 40-bit accumulator, shifted right by shift
        and saturated to 32 bits.
        \param[in] in1 vector (16 bytes aligned, 16 readable bytes after the end).
        \param[in] in2 vector (16 bytes aligned, 16 readable bytes after the end).
        \param[in] size vector size (multiple of 8, >= 8).
        \param[in] shift right shift of the sum (< 40).
        \return (sum in1[i] * in2[i]) >> shift, saturated to int32
    */
    int32_t dot_product_16_16_32(q15 *in1, q15 *in2, uint32_t size, uint32_t shift);
    /// Dot product of vectors q15.
    /*!
        \param[in] in1 vector.
        \param[in] in2 vector (16 bytes aligned).
        \param[in] size vector size (multiple of 8, >= 8).
        \return result
    */
    q15 dot_product_1_16(q15 *in1, q15 *in2, uint32_t size);
    /// FIR q15.
    /*!
        \param[in] in data.
        \param[in] k coefficients (16 bytes aligned).
        \param[in] ksize size of coefficients (multiple of 8, >= 8).
        \param[out] out output vector.
        \param[in] size vector size.
    */
    void fir_1_16(q15 *in, q15 *k, uint32_t ksize, q15 *out, uint32_t size);
    /// Sums of sample pairs on 12 adjacent shifts (sparse correlation with +-1 taps).
    /*!
        out[k] = sum over j of (in[pos[2j+1] + k] - in[pos[2j] + k]), k = 0..11, exact (32 bit).
        \param[in] in samples (any alignment); in[p + 0..23] must be readable for every position p.
        \param[in] pos positions of the samples: minus, plus, minus, plus, ...
        \param[in] pairs number of pairs.
        \param[out] out 12 sums (16 bytes aligned).
    */
    void sum_pairs_12_s16(const int16_t *in, const uint16_t *pos, uint32_t pairs, int32_t *out);

    /// Addition of two vectors with saturation.
    /*!
        out=in1 + in2
        \param[in] in1 q15 vector (16 bytes aligned).
        \param[in] in2 q15 vector (16 bytes aligned).
        \param[out] out sum vector (16 bytes aligned).
        \param[in] size vector size (multiple of 8).
    */
    void addVectors_q15(q15 *in1, q15 *in2, q15 *out, uint32_t size);

    /// Normalize vector q15.
    /*!
        \param[in] in q15 vector (16 bytes aligned).
        \param[out] out q15 vector (16 bytes aligned).
        \param[in] size vector size (multiple of 8).
        \return 0,1 if normalization is not needed
    */
    int16_t normalize_q15(q15 *in, q15 *out, uint32_t size);

    /// Normalize vector q14.
    /*!
        \param[in] in q15 vector (16 bytes aligned).
        \param[out] out q15 vector (16 bytes aligned).
        \param[in] size vector size (multiple of 8).
        \return 0,2,3 if normalization is not needed
    */
    int16_t normalize_q14(q15 *in, q15 *out, uint32_t size);

    /// FIR q15.
    /*!
        \param[in] in data vector (16 bytes aligned).
        \param[in] k coefficients (16 bytes aligned).
        \param[in] ksize size of coefficients (multiple of 8, >= 8).
        \param[out] out output vector (16 bytes aligned).
        \param[in] size vector size (multiple of 8, >= 8).
    */
    void fir_16_16_q15(q15 *in, q15 *k, uint32_t ksize, q15 *out, uint32_t size);

    /// atan2 q15.
    /*!
        \param[in] y.
        \param[in] x.
        \return atan(y/x) in radians. (Pi = 16383)
    */
    q15 atan2_q15(q15 y, q15 x);

    /// atan q15.
    /*!
        \param[in] y tangent.
        \return atan(y) in radians. (Pi = 16383)
    */
    inline q15 atan_q15(q15 y)
    {
        return atan2_q15(y, 0x7fff);
    };

    /// sin cos q15.
    /*!
        \param[in] angle angle. (Pi = 16383)
        \param[out] sn sine.
        \param[out] cs cosine.
    */
    void sincos_q15(q15 angle, q15 *sn, q15 *cs);

#ifdef __cplusplus
}
#endif