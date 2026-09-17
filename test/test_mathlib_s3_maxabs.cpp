/**
 * @file test_mathlib_s3_maxabs.cpp
 * @author Roman Bliznets (r.bliznets@gmail.com)
 * @brief Unity module tests for maxAbsVector_16 (IQ 24->16 shift calculation)
 * @version 0.0.0.1
 * @date 17.09.2026
 *
 * @copyright Copyright (c) 2026
 *
 */

#include <limits.h>
#include <cstring>
#include "unity.h"
#include "mathlib_s3.h"
#include "CTrace.h"
#include <math.h>
#include "esp_system.h"
#include "esp_random.h"

#define countof(x) (sizeof(x)/sizeof(x[0]))
#define N       10000   /**< Number of iterations for benchmark */
#define WORDS   24      /**< Size of the packet in words: to8SIZE(RX_SYM_NUM*2) = to8SIZE(20) = 24 */
#define DATA_WORDS 20   /**< Number of words written by readIQ: RX_SYM_NUM*2 = 20 */

__attribute__((aligned(16)))
static uint32_t packet[WORDS];

/**
 * @brief Reference scalar implementation of the shift calculation (former CI2STask code)
 *
 * @param data raw IQ words in AK2401 format (bit 31 - sign, bits 30..7 - 24bit sample)
 * @return shift value
 */
static uint8_t ref_shift(const uint32_t* data)
{
    uint8_t shift = 7;
    int32_t mx = 0;
    for (int i = 0; i < DATA_WORDS; i++)
    {
        int32_t v = (int32_t)(data[i] << 1); ///< Shift by 1 bit (data format).
        if (v < 0)
            v = (v == INT32_MIN) ? INT32_MAX : -v; ///< Absolute value (INT32_MIN saturates).
        if (v > mx) mx = v;
    }
    mx >>= 7; ///< Normalize the maximum relative to the initial shift.
    while (mx > INT16_MAX)
    {
        mx >>= 1;
        shift++;
    }
    return shift;
}

/**
 * @brief SIMD implementation of the shift calculation (new CI2STask code)
 *
 * @param data raw IQ words in AK2401 format
 * @return shift value
 */
static uint8_t simd_shift(const uint32_t* data)
{
    uint8_t shift = 7;
    uint32_t mx = maxAbsVector_16((uint32_t*)data, WORDS);
    while (mx >= 32)
    {
        mx >>= 1;
        shift++;
    }
    return shift;
}

/**
 * @brief Reference scalar implementation of maxAbsVector_16 (C equivalent)
 *
 * @param data input words
 * @param size number of words
 * @return max abs of (int16_t)(data[i] >> 16)
 */
static uint32_t ref_maxAbs(const uint32_t* data, uint32_t size)
{
    int16_t mn = 0, mx = 0;
    for (uint32_t i = 0; i < size; i++)
    {
        int16_t v = (int16_t)(((int32_t)data[i]) >> 16); ///< High half of the word (arithmetic shift).
        if (v < mn) mn = v;
        if (v > mx) mx = v;
    }
    return (uint32_t)((mx >= -mn) ? mx : -mn);
}

TEST_CASE("maxAbsVector_16", "[math][mathlib_s3]")
{
    uint32_t mem1 = esp_get_free_heap_size();

    TRACE("Number of iterations", N, false);
    TRACE("Packet size (words)", WORDS, false);

    /// Direct check of H against the C equivalent on all boundary patterns.
    /// Patterns: zeros, positive/negative full-scale sample, asymmetric extremes,
    /// -32768 in the high half (min reduction must survive abs), random noise.
    const uint32_t patterns[][4] = {
        {0, 0, 0, 0},                                        // silence
        {0x7FFFFF80u, 0x7FFFFF80u, 0x7FFFFF80u, 0x7FFFFF80u},// positive full-scale sample
        {0x80000080u, 0x80000080u, 0x80000080u, 0x80000080u},// negative full-scale sample
        {0xFF7FFF80u, 0x00800080u, 0x7F000080u, 0x81000080u},// -32768/+32768 in high halves
        {0x00008000u, 0xFFFF8000u, 0x12345678u, 0xFEDCBA98u},// high half exactly 0x8000 (=> -32768)
    };
    for (auto& p : patterns)
    {
        for (int i = 0; i < WORDS; i++)
            packet[i] = p[i & 3];
        uint32_t expect = ref_maxAbs(packet, WORDS);
        uint32_t got = maxAbsVector_16(packet, WORDS);
        if (got != expect)
        {
            TDEC("expected H", expect);
            TDEC("got H", got);
            TEST_FAIL_MESSAGE("maxAbsVector_16 != C reference");
        }
    }

    /// Random packets in AK2401 format: w = (uint32_t)((int32_t)s << 7) | (status & 0x7f),
    /// s is a 24bit sample covering all amplitudes (including full-scale ±2^23).
    /// Checked: simd shift is never less than reference and at most +1 (safe side margin),
    /// and the output word after >> shift fits into int16_t without wraparound.
    for (int t = 0; t < 2000; t++)
    {
        int32_t amp = 1 << (1 + (t % 23)); ///< Amplitude sweeps all bit widths 2^1..2^23.
        for (int i = 0; i < DATA_WORDS; i++)
        {
            int32_t s = (int32_t)(((esp_random() >> 8) % (2 * amp + 1)) - amp); ///< Uniform in [-amp, amp].
            packet[i] = (uint32_t)((s << 7)) | (esp_random() & 0x7f);
        }
        for (int i = DATA_WORDS; i < WORDS; i++)
            packet[i] = 0; ///< Padding beyond 20 words must be zeroed (as in CI2STask MSG_START_IQ).

        uint8_t ref = ref_shift(packet);
        uint8_t sim = simd_shift(packet);
        if ((sim < ref) || (sim > ref + 1))
        {
            TDEC("test", t);
            TDEC("ref shift", ref);
            TDEC("simd shift", sim);
            TEST_FAIL_MESSAGE("simd shift out of {ref, ref+1}");
        }

        /// No int16 wraparound after the selected shift: |(int32)(w>>1)>>sim| <= INT16_MAX.
        for (int i = 0; i < DATA_WORDS; i++)
        {
            int32_t v = (int32_t)(packet[i] << 1) >> sim;
            if ((v > INT16_MAX) || (v < -INT16_MAX))
            {
                TDEC("test", t);
                TDEC("shift", sim);
                TEST_FAIL_MESSAGE("int16 overflow after simd shift");
            }
        }
    }

    /// Benchmark: scalar packet scan (former CI2STask loop) vs PIE.
    for (int i = 0; i < DATA_WORDS; i++)
        packet[i] = (uint32_t)((int32_t)(i * 700000) << 7) | (i & 0x7f);
    for (int i = DATA_WORDS; i < WORDS; i++)
        packet[i] = 0;

    /// The volatile accumulator plus the compiler barrier keep both loops alive:
    /// ref_shift is a side-effect-free static function, so without an observable
    /// use of its result GCC deletes the scalar loop (it measured 0 nsec).
    static volatile uint8_t sink = 0;
    STARTTIMESHOT();
    for (int i = 0; i < N; i++)
    {
        sink += ref_shift(packet);
        asm volatile("" ::: "memory");
    }
    STOPTIME("scalar shift calc time", N);

    STARTTIMESHOT();
    for (int i = 0; i < N; i++)
    {
        sink += simd_shift(packet);
        asm volatile("" ::: "memory");
    }
    STOPTIME("maxAbsVector_16 shift calc time", N);

    uint32_t mem2 = esp_get_free_heap_size();
    if (mem1 != mem2)
    {
        TRACE("memory leak", mem1 - mem2, false);
        TRACE("start", mem1, false);
        TRACE("stop", mem2, false);
        TEST_FAIL_MESSAGE("memory leak");
    }
}
