# Q15 Algorithm Library Optimized for esp32-s3
To add to a project in the components folder from the command line, run:    

    git submodule add https://github.com/rbliznets/esp32-mathlib_s3 mathlib_s3
## FFT
Radix2 implementation with and without auto-normalization. Execution time:
```
This is esp32s3 (240MHz) chip with 2 CPU core(s), WiFi/BLE, silicon revision 1, Minimum free heap size: 378624 bytes

**** FFT size 8192 *****
fft_radix2(1): 484usec
fft_radix2_scale(1): 904usec
revbin_permute(1): 260usec

**** FFT size 4096 *****
fft_radix2(1): 222usec
fft_radix2_scale(1): 376usec
revbin_permute(1): 130usec

**** FFT size 2048 *****
fft_radix2(1): 107usec
fft_radix2_scale(1): 183usec
revbin_permute(1): 66usec
```

Library contains optimized functions for performing Fast Fourier Transforms (FFT) on the ESP32-S3, operating on fixed-point Q15 complex numbers.

*   **`fft_log2(uint32_t fftSize)`**:
    *   Calculates the base-2 logarithm of the given FFT size.
    *   Returns the calculated log2 value.

*   **`init_fft(complex_q15 *w, uint32_t fftSize)`**:
    *   Initializes a lookup table (`w`) containing the complex twiddle factors required for the FFT computation.
    *   The table is sized for the maximum `fftSize` specified.

*   **`getW(complex_q15 *w, uint32_t fftSize, uint32_t fftSize2)`**:
    *   Retrieves a pointer to the section of the pre-calculated twiddle factor table (`w`) that corresponds to a specific `fftSize`.
    *   The full table `w` was initialized for a maximum size of `fftSize2`.

*   **`revbin_permute(complex_q15 *data, uint32_t fftSize)`**:
    *   Performs an in-place bit-reversal permutation on the `data` array.
    *   This step is typically required after an FFT algorithm to reorder the output into the correct sequence.

*   **`fft_radix2(complex_q15 *data, complex_q15 *w, uint32_t fftSize)`**:
    *   Performs an in-place Radix-2 FFT on the `data` array.
    *   Uses the twiddle factors from array `w`.
    *   Applies a fixed scaling factor of 1/2 at each stage of the FFT computation.

*   **`fft_radix2_scale(complex_q15 *data, complex_q15 *w, uint32_t fftSize)`**:
    *   Performs an in-place Radix-2 FFT on the `data` array.
    *   Uses the twiddle factors from array `w`.
    *   Applies automatic scaling during the computation to help prevent overflow.

## CORDIC
Library contains defines optimized functions for performing operations on fixed-point Q15 complex numbers, specifically designed for the ESP32-S3.

*   **`atan2_q15(q15 y, q15 x)`**:
    *   Calculates the 4-quadrant arc tangent of `y/x`.
    *   Returns the angle in radians, scaled such that π corresponds to 16383.

*   **`atan_q15(q15 y)`**:
    *   An inline function that calculates the arc tangent of `y` (tangent value).
    *   It uses `atan2_q15(y, 0x7fff)` internally, effectively treating `y` as `y/1`.
    *   Returns the angle in radians, scaled such that π corresponds to 16383.

*   **`sincos_q15(q15 angle, q15 *sn, q15 *cs)`**:
    *   Calculates the sine and cosine of a given `angle`.
    *   The `angle` is provided in the scaled Q15 format where π corresponds to 16383.
    *   The resulting sine value is stored in `sn`, and the cosine value is stored in `cs`.

*   **`arg(complex_q15 value)`**:
    *   An inline function that calculates the argument (angle/phase) of a single `complex_q15` value.
    *   It uses the `atan2_q15` function internally.
    *   Returns the argument in radians, scaled such that π corresponds to 16383.

*   **`arg_16_q15(complex_q15 *in, q15 *out, uint32_t size)`**:
    *   Calculates the argument (angle/phase) for each complex number in an input vector (`in`).
    *   The input vector must be 16-byte aligned, and its `size` must be a multiple of 8.
    *   The resulting arguments are stored in the output vector (`out`), which must also be 16-byte aligned.

Execution time:
```
sincos_q15(1000): 915nsec
**** Data size 1016 *****
arg_q15(1): 214usec
arg_q15(1000): 198usec
**** Data size 1 *****
float std::arg(1000): 1096nsec
arg_fr16(1000): 13usec
arg(1000): 709nsec
```
## Complex
Library contains defines optimized functions for performing operations on fixed-point Q15 complex numbers, specifically designed for the ESP32-S3.

*   **`complex_q15`**:
    *   A structure representing a complex number using two `q15` (16-bit fixed-point) values for the real (`re`) and imaginary (`im`) parts. It is aligned to 4 bytes.

*   **`magnitude_q15(complex_q15 *in, q15 *out, uint32_t size)`**:
    *   Calculates the square magnitude divided by 2 (`(real^2 + imag^2) / 2`) for each complex number in an input vector (`in`).
    *   The input vector must be 16-byte aligned, and its `size` must be a multiple of 8.
    *   The resulting magnitudes are stored in the output vector (`out`), which must also be 16-byte aligned.

*   **`cmul_q15(complex_q15 x, complex_q15 y)`**:
    *   Performs complex multiplication of two `complex_q15` numbers (`x` and `y`).
    *   Returns the resulting complex number.

*   **`cmul10_q15(complex_q15 *in, complex_q15 *k, complex_q15 *out)`**:
    *   Multiplies each element of a complex input vector (`in`, implicitly of size 10) by a single complex scalar (`k`).
    *   The input vector and output vector (`out`) must be 16-byte aligned.
    *   The scalar `k` must be 2-byte aligned.
    *   Stores the results in the output vector (`out`).

## Functions
*   **`copy(q15 *in, q15 *out, uint32_t size)`**:
    *   Copies `size` elements from the input vector `in` to the output vector `out`.
    *   The output vector `out` must be 16-byte aligned. The `size` must be a multiple of 8 and >= 8.

*   **`copy_16(q15 *in, q15 *out, uint32_t size)`**:
    *   Copies `size` elements from the input vector `in` to the output vector `out`.
    *   Both input `in` and output `out` vectors must be 16-byte aligned. The `size` must be a multiple of 8 and >= 8.

*   **`scaleVector_16(q15 *in, q15 *k, q15 *out, uint32_t size)`**:
    *   Multiplies each element of the input vector `in` by a scalar value pointed to by `k`, storing the result in `out`.
    *   The input `in` and output `out` vectors must be 16-byte aligned. The scalar pointer `k` must be 2-byte aligned. The `size` must be a multiple of 8 and >= 8.

*   **`shrinkVector_16(uint32_t *in, uint8_t shift, q15 *out, uint32_t size)`**:
    *   Performs an arithmetic right shift (`shift` bits) on each element of a 32-bit input vector `in` and stores the upper 15 bits of the result in the 16-bit output vector `out`.
    *   The input `in` and output `out` vectors must be 16-byte aligned. The `shift` value must be between 1 and 31. The `size` must be a multiple of 8 and >= 8.

*   **`maxAbsVector_16(uint32_t *in, uint32_t size)`**:
    *   Finds the maximum absolute value of `(int16_t)(in[i] >> 16)` across all elements of the 32-bit input vector `in` (i.e. the high halves of the words, as produced by the AK2401 IQ format).
    *   The input `in` vector must be 16-byte aligned. The `size` (in words) must be a multiple of 8 and >= 8.
    *   Returns the maximum absolute value (0..32768).

*   **`scaleVector(q15 *in, q15 *k, q15 *out, uint32_t size)`**:
    *   Multiplies each element of the input vector `in` by a scalar value pointed to by `k`, storing the result in `out`.
    *   The scalar pointer `k` must be 2-byte aligned, and the output `out` must be 16-byte aligned. The `size` must be a multiple of 8 and >= 16.

*   **`addVectors_q15(q15 *in1, q15 *in2, q15 *out, uint32_t size)`**:
    *   Performs element-wise addition of two vectors `in1` and `in2`, storing the result in `out` with saturation to prevent overflow.
    *   All vectors (`in1`, `in2`, `out`) must be 16-byte aligned. The `size` must be a multiple of 8.

*   **`normalize_q15(q15 *in, q15 *out, uint32_t size)`**:
    *   Normalizes the input vector `in` such that the maximum absolute value fits within the Q15 range and stores the result in `out`.
    *   Both vectors must be 16-byte aligned. The `size` must be a multiple of 8.
    *   Returns 0 or 1 if normalization was not necessary (based on the maximum value found), otherwise returns a value indicating normalization was performed.

*   **`normalize_q14(q15 *in, q15 *out, uint32_t size)`**:
    *   Normalizes the input vector `in` such that the maximum absolute value fits within the Q14 range (using an internal Q15 representation) and stores the result in `out`.
    *   Both vectors must be 16-byte aligned. The `size` must be a multiple of 8.
    *   Returns 0, 2, or 3 if normalization was not necessary, otherwise returns a value indicating normalization was performed.

## Filter
Library contains defines optimized functions for performing operations on fixed-point Q15 complex numbers, specifically designed for the ESP32-S3.

*   **`dot_product_16_16(q15 *in1, q15 *in2, uint32_t size)`**:
    *   Calculates the dot product (sum of element-wise products) of two vectors `in1` and `in2`.
    *   Both input vectors must be 16-byte aligned. The `size` must be a multiple of 8 and >= 8.
    *   Returns the resulting scalar value.

*   **`dot_product_1_16(q15 *in1, q15 *in2, uint32_t size)`**:
    *   Calculates the dot product of two vectors `in1` and `in2`.
    *   The second input vector `in2` must be 16-byte aligned. The `size` must be a multiple of 8 and >= 8. The first vector `in1` has no specified alignment requirement in this signature's comment.
    *   Returns the resulting scalar value.
  
*   **`fir_1_16(q15 *in, q15 *k, uint32_t ksize, q15 *out, uint32_t size)`**:
    *   Performs Finite Impulse Response (FIR) filtering on input data `in` using coefficients `k`.
    *   The coefficients `k` must be 16-byte aligned. The number of coefficients `ksize` must be a multiple of 8 and >= 8.
    *   Outputs `size` samples to the `out` vector.

*   **`fir_16_16_q15(q15 *in, q15 *k, uint32_t ksize, q15 *out, uint32_t size)`**:
    *   Performs Finite Impulse Response (FIR) filtering on input data `in` using coefficients `k`.
    *   All vectors (`in`, `k`, `out`) must be 16-byte aligned. The number of coefficients `ksize` must be a multiple of 8 and >= 8. The output `size` must be a multiple of 8 and >= 8.
    *   Outputs `size` samples to the `out` vector.

Execution time:
```
**** Data size 1016 *****
dot_product_16_16(1000): 1765nsec
dot_product_1_16(1000): 2310nsec
dot_product_c(1000): 38usec

**** Data size 1000 *****
fir_1_16(1000): 196usec
fir_16_16_q15(1000): 36usec
```

## Definitions and Macros:

*   **`q15`**: A typedef for `int16_t`, representing a signed fixed-point number with 1 bit for sign and 15 bits for magnitude (S1.15 format).
*   **`toQ15(f)`**: Macro to convert a float `f` to `q15` format.
*   **`toFloat(q)`**: Macro to convert a `q15` value `q` back to float.
*   **`toQ15Angle(f)`**: Macro to convert a float angle `f` (in radians) to `q15` format, scaled so that 2π radians maps to the full range of `q15`.
*   **`toFloatAngle(q)`**: Macro to convert a `q15` angle `q` back to float radians.
*   **`to16SIZE(size)`**: Macro to round up `size` to the nearest multiple of 16.
*   **`to8SIZE(size)`**: Macro to round up `size` to the nearest multiple of 8.
*   **`to4SIZE(size)`**: Macro to round up `size` to the nearest multiple of 4.

