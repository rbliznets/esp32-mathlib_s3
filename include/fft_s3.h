/*!
	\file
	\brief Оптимизированные функции для FFT.
	\authors Близнец Р.А.
	\version 0.0.0.1
	\date 13.03.2023
*/

#pragma once
#include "complex_s3.h"

#ifdef __cplusplus
extern "C"
{
#endif

/// log2 for FFT.
/*!
    \param[in] fftSize FFT size.
    \return log2
*/
int16_t fft_log2(uint32_t fftSize);
/// Init FFT twidle.
/*!
	Initialization FFT twidle.
    \param[out] w twiddle array (fftSize)
    \param[in] fftSize maximum FFT size.
*/
void init_fft(complex_q15* w, uint32_t fftSize);
/// get FFT twidle array.
/*!
    \param[in] w twiddle array (fftSize2)
    \param[in] fftSize FFT size.
    \param[in] fftSize2 maximum FFT size.
    \return begin FFT twidle array for fftSize.
*/
complex_q15* getW(complex_q15* w, uint32_t fftSize, uint32_t fftSize2);

/// permute data after FFT.
/*!
    \param[in|out] data FFT output
    \param[in] fftSize FFT size.
*/
void revbin_permute(complex_q15* data, uint32_t fftSize);

/// FFT.
/*!
	Stage scailing: 1/2
    \param[in|out] data 
    \param[in] w twiddle array (fftSize)
    \param[in] fftSize FFT size.
*/
void fft_radix2(complex_q15* data, complex_q15* w, uint32_t fftSize);
/// FFT.
/*!
	Stage scailing: auto
    \param[in|out] data 
    \param[in] w twiddle array (fftSize)
    \param[in] fftSize FFT size.
*/
void fft_radix2_scale(complex_q15* data, complex_q15* w, uint32_t fftSize);

#ifdef __cplusplus
}
#endif

