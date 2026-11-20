/*!
    \file
    \brief Optimized functions for complex numbers.
    \authors Bliznets R.A. (r.bliznets@gmail.com)
    \version 1.0.0.0
    \date 02/14/2023
*/

#pragma once
#include "mathlib_s3.h"

#ifdef __cplusplus
extern "C"
{
#endif
    /// Complex number.
    typedef __attribute__((aligned(4))) struct complex_q15
    {
        q15 re; ///< real
        q15 im; ///< imaginary
    } complex_q15;

    /// argument.
    /*!
        \param[in] value The complex number.
        \return argument in radians. (Pi = 16383)
    */
    inline q15 arg(complex_q15 value)
    {
        return atan2_q15(value.im, value.re);
    };

    /// argument.
    /*!
        \param[in] in complex vector (16 bytes aligned).
        \param[out] out argument in radians (16 bytes aligned).
        \param[in] size vector size (multiple of 8).
    */
    void arg_16_q15(complex_q15 *in, q15 *out, uint32_t size);

    /// The square magnitude of a complex vector divided by 2.
    /*!
        out=(in.re^2 + in.im^2)/2
        \param[in] in complex vector (16 bytes aligned).
        \param[out] out magnitude (16 bytes aligned).
        \param[in] size vector size (multiple of 8).
    */
    void magnitude_q15(complex_q15 *in, q15 *out, uint32_t size);

    /// Complex multiplication.
    /*!
        out = x * y
        \param[in] x complex value.
        \param[in] y complex value.
        \return x*y.
    */
    complex_q15 cmul_q15(complex_q15 x, complex_q15 y);
    /// Multiply complex vector (size=10) by complex scalar.
    /*!
        out=in * k
        \param[in] in vector (16 bytes aligned).
        \param[in] k pointer to scalar (2 bytes aligned).
        \param[out] out output vector (16 bytes aligned).
    */
    void cmul10_q15(complex_q15 *in, complex_q15 *k, complex_q15 *out);
    /// Rotate a complex vector by a periodic twiddle table.
    /*!
        data[i] = data[i] * tw[(phase + i) % period]
        The table is stored per phase: tw[p][j] is the twiddle of the phase (p + j) % period, so four
        consecutive twiddles are always read with one 128-bit load even if period is not a multiple of 4.
        \param[in,out] data vector (16 bytes aligned), modified in place.
        \param[in] size vector size (multiple of 4).
        \param[in] tw phase table (16 bytes aligned).
        \param[in] period number of phases (>= 1).
        \param[in] phase starting phase (0..period-1).
    */
    void rotate_q15(complex_q15 *data, uint32_t size, const complex_q15 (*tw)[4], uint32_t period, uint32_t phase);

    /// Split an interleaved complex vector into real and imaginary parts.
    /*!
        re[i] = in[i].re, im[i] = in[i].im
        \param[in] in complex vector (16 bytes aligned).
        \param[out] re real part (16 bytes aligned).
        \param[out] im imaginary part (16 bytes aligned).
        \param[in] size vector size (multiple of 8).
    */
    void iq_unzip_q15(complex_q15 *in, q15 *re, q15 *im, uint32_t size);

    /// Merge real and imaginary parts into an interleaved complex vector.
    /*!
        out[i].re = re[i], out[i].im = im[i]
        \param[in] re real part (16 bytes aligned).
        \param[in] im imaginary part (16 bytes aligned).
        \param[out] out complex vector (16 bytes aligned).
        \param[in] size vector size (multiple of 8).
    */
    void iq_zip_q15(q15 *re, q15 *im, complex_q15 *out, uint32_t size);

#ifdef __cplusplus
}
#endif