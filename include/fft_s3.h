/*!
    \file
    \brief Optimized functions for FFT.
    \authors Bliznets R.A. (r.bliznets@gmail.com)
    \version 1.0.0.0
    \date 03/13/2023
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
    /// Initialize FFT twiddle factors.
    /*!
        Initializes the FFT twiddle factor array.
        \param[out] w twiddle array (fftSize elements)
        \param[in] fftSize maximum FFT size.
    */
    void init_fft(complex_q15 *w, uint32_t fftSize);
    /// Get FFT twiddle array for a specific size.
    /*!
        \param[in] w twiddle array (fftSize2 elements)
        \param[in] fftSize desired FFT size.
        \param[in] fftSize2 maximum FFT size (size of the w array).
        \return pointer to the start of the FFT twiddle array for fftSize.
    */
    complex_q15 *getW(complex_q15 *w, uint32_t fftSize, uint32_t fftSize2);

    /// Permute data after FFT (bit-reversal).
    /*!
        \param[in|out] data FFT output
        \param[in] fftSize FFT size.
    */
    void revbin_permute(complex_q15 *data, uint32_t fftSize);

    /// FFT.
    /*!
        Stage scaling: 1/2
        \param[in|out] data
        \param[in] w twiddle array (fftSize elements)
        \param[in] fftSize FFT size.
    */
    void fft_radix2(complex_q15 *data, complex_q15 *w, uint32_t fftSize);
    /// FFT.
    /*!
        Stage scaling: automatic
        \param[in|out] data
        \param[in] w twiddle array (fftSize elements)
        \param[in] fftSize FFT size.
    */
    void fft_radix2_scale(complex_q15 *data, complex_q15 *w, uint32_t fftSize);

#ifdef __cplusplus
}
#endif