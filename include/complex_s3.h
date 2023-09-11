/*!
	\file
	\brief Оптимизированные функции для комплексных чисел.
	\authors Близнец Р.А.
	\version 0.3.0.0
	\date 14.02.2023
*/

#pragma once
#include "mathlib_s3.h"

#ifdef __cplusplus
extern "C"
{
#endif
/// Complex number.
typedef __attribute__((aligned(4))) struct complex_q15{
    q15 re; ///< real
    q15 im; ///< imagine
} complex_q15;

/// argument.
/*!
    \param[in] value .
    \return argument in radians. (Pi = 16383)
*/
inline q15 arg(complex_q15 value)
{
    return atan2_q15(value.im, value.re);
};

/// argument.
/*!
    \param[in] complex vector (16 bytes aligned)..
    \param[out] out magnitude (16 bytes aligned).
    \param[in] size vector size (multiple of 8).
*/
void arg_16_q15(complex_q15* in, q15* out, uint32_t size);

/// The square magnitude of a complex vector divided by 2.
/*!
    out=(in.re^2 + in.im^2)/2
    \param[in] in complex vector (16 bytes aligned).
    \param[out] out magnitude (16 bytes aligned).
    \param[in] size vector size (multiple of 8).
*/
void magnitude_q15(complex_q15* in, q15* out, uint32_t size);

/// Complex multiplication.
/*!
    out=(in.re^2 + in.im^2)/2
    \param[in] x complex value.
    \param[out] y complex value.
    \return x*y.
*/
complex_q15 cmul_q15(complex_q15 x, complex_q15 y);
/// Multiply complex vector (size=10) by complex scalar.
/*!
    out=in * k
    \param[in] in vector (16 bytes aligned).
    \param[in] k pointer for scalar (2 bytes aligned).
    \param[out] out output vector (16 bytes aligned).
*/
void cmul10_q15(complex_q15* in, complex_q15* k, complex_q15* out);

#ifdef __cplusplus
}
#endif
