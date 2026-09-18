/*!
    \file
    \brief Optimized functions for FFT.
    \authors Bliznets R.A. (r.bliznets@gmail.com)
    \version 1.0.0.0
    \date 03/13/2023
*/

#include "fft_s3.h"
#include <math.h>

#include <stdio.h>
#include "esp_attr.h"

// Calculates the base-2 logarithm of the FFT size.
// Returns the log2 value if fftSize is a power of two, otherwise returns a negative value.
int16_t fft_log2(uint32_t fftSize)
{
    uint32_t x = 0x8000; // Start with the highest bit for a 16-bit value
    int16_t res = 15;    // Corresponding log2 value
    while (x != 0)
    {
        if (((x & fftSize) != 0) || (res == 0))
        {
            if (x == fftSize) // fftSize is a power of two
            {
                return res;
            }
            else // fftSize is not a power of two
            {
                return -res;
            }
        }
        res--;   // Decrement the log2 estimate
        x >>= 1; // Shift the bit mask right
    }
    return 0; // Should not be reached for valid inputs
}

// Updates the bit-reversed index for the next iteration in bit-reversal permutation.
inline uint32_t revbin_update(uint32_t r, uint32_t n)
{
    // This loop manipulates 'r' to get the next bit-reversed index based on 'n'
    for (uint32_t m = n >> 1; (!((r ^= m) & m)); m >>= 1)
        ;
    return r;
}

#if defined(__XTENSA__)
// Internal implementation of bit-reversal permutation for sizes 16..1024 (PIE, EE.BITREV)
void revbin_permute_pie(complex_q15 *data, uint32_t fftSize, uint32_t width);
#endif

// Performs bit-reversal permutation on the FFT data array in-place.
void IRAM_ATTR revbin_permute(complex_q15 *data, uint32_t fftSize)
{
#if defined(__XTENSA__)
    // EE.BITREV covers 4..10 bits (sizes 16..1024)
    if ((fftSize >= 16) && (fftSize <= 1024))
    {
        revbin_permute_pie(data, fftSize, ((uint32_t)fft_log2(fftSize)) & 7);
        return;
    }
#endif
    uint32_t *dt = (uint32_t *)data; // Treat complex_q15 pairs as 32-bit words for faster swapping
    uint32_t nh = fftSize >> 1;      // Half the FFT size
    uint32_t r = 0;                  // Bit-reversed index
    uint32_t x = 1;                  // Forward index
    uint32_t t;                      // Temporary variable for swapping
    while (x < nh)
    {
        r = r + nh; // Update r for the next pair of swaps
        t = dt[x];  // Swap elements at positions x and r
        dt[x] = dt[r];
        dt[r] = t;
        x++; // Move to the next pair

        r = revbin_update(r, fftSize); // Calculate next bit-reversed index
        if (r > x)                     // Only swap if r is greater than x to avoid double-swapping
        {
            t = dt[x]; // Swap elements at positions x and r
            dt[x] = dt[r];
            dt[r] = t;
            t = dt[fftSize - 1 - x]; // Swap corresponding elements from the end
            dt[fftSize - 1 - x] = dt[fftSize - 1 - r];
            dt[fftSize - 1 - r] = t;
        }
        x++; // Move to the next pair
    }
}

// Initializes the twiddle factor lookup table for the FFT.
// fftSize - 2 (This comment seems incorrect, it should likely be fftSize/2 or similar)
void init_fft(complex_q15 *w, uint32_t fftSize)
{
    int16_t n = fft_log2(fftSize);

    assert(((uint32_t)w % 16) == 0); // Ensure alignment
    assert(n >= 3);                  // Ensure minimum size

    float e = M_PI * 2.0 / fftSize; // Fundamental angle for twiddle factors
    fftSize >>= 1;                  // Calculate for half the size initially
    for (int i = 0; i < fftSize; i++)
    {
        // Calculate cosine and sine components and convert to Q15 format
        w[i].re = (q15)roundf(INT16_MAX * cosf(i * e));
        w[i].im = (q15)roundf(INT16_MAX * sinf(i * e));
    }

    // Generate twiddle factors for smaller stages by subsampling the largest stage
    complex_q15 *w_last = w;               // Pointer to the previous stage's factors
    complex_q15 *w_cur = &w_last[fftSize]; // Pointer for the current stage
    for (int i = 3; i < n; i++)
    {
        fftSize >>= 1; // Halve the size for the next stage
        for (int j = 0; j < fftSize; j++)
        {
            // Copy every second factor from the previous stage
            w_cur[j].re = w_last[2 * j].re;
            w_cur[j].im = w_last[2 * j].im;
        }
        w_last = w_cur; // Move pointers to the current stage for the next iteration
        w_cur = &w_last[fftSize];
    }
    // Handle the final stage (size 2) by negating the factors from the previous stage
    fftSize >>= 1;
    for (int j = 0; j < fftSize; j++)
    {
        w_cur[j].re = -w_last[2 * j].re;
        w_cur[j].im = -w_last[2 * j].im;
    }
}

// Gets the pointer to the twiddle factor array for a specific FFT size from the pre-calculated table.
complex_q15 *getW(complex_q15 *w, uint32_t fftSize, uint32_t fftSize2)
{
    assert(fftSize >= fftSize2); // fftSize2 should not be larger than the maximum fftSize used in init_fft

    complex_q15 *w_cur = w; // Start at the beginning of the table
    while (fftSize > fftSize2)
    {
        fftSize >>= 1;           // Navigate down the table structure
        w_cur = &w_cur[fftSize]; // Move to the next level
    }
    return w_cur; // Return the pointer to the requested size's factors
}

// Internal implementation function for radix-2 FFT (PIE - Peripheral Input/Output Engine or specific hardware acceleration)
void fft_r2_q15_pie(complex_q15 *data, complex_q15 *w, uint32_t fftSize);
// Wrapper function for radix-2 FFT with stage scaling 1/2, including assertions
inline void fft_radix2(complex_q15 *data, complex_q15 *w, uint32_t fftSize)
{
    assert(((uint32_t)data % 16) == 0); // Ensure alignment
    assert(fftSize >= 16);              // Ensure minimum size

    fft_r2_q15_pie(data, w, fftSize);
}

// Internal implementation function for scaled radix-2 FFT (PIE)
void fft_r2s_q15_pie(complex_q15 *data, complex_q15 *w, uint32_t fftSize);
// Wrapper function for scaled radix-2 FFT with automatic scaling, including assertions
inline void fft_radix2_scale(complex_q15 *data, complex_q15 *w, uint32_t fftSize)
{
    assert(((uint32_t)data % 16) == 0); // Ensure alignment
    assert(fftSize >= 16);              // Ensure minimum size

    fft_r2s_q15_pie(data, w, fftSize);
}