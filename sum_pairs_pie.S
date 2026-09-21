    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access)
    .align  4
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Sums of sample pairs at arbitrary positions on 12 adjacent shifts (sparse correlation with +-1 taps):
//     out[k] = sum over j < pairs of (in[pos[2j+1] + k] - in[pos[2j] + k]), k = 0..11.
//
// For each position 16 samples are read: three aligned 128-bit loads (EE.LD.128.USAR.IP sets SAR_BYTE
// from the low address bits, the same for all three) and two funnel shifts EE.SRC.Q. The samples are
// sign-extended to 32 bits (EE.VCMP.LT.S16 against zero gives the high halves, EE.VZIP.16 interleaves
// them) and accumulated with EE.VADDS.S32 / EE.VSUBS.S32, so the sums are exact while they fit in
// 32 bits (at most 2^16 pairs of full-scale samples).
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   sum_pairs_12_s16_pie,@function
    .align  4
    .global sum_pairs_12_s16_pie
sum_pairs_12_s16_pie:
// in    - a2  (pointer to int16_t samples, any alignment; in[pos + 0..23] must be readable)
// pos   - a3  (pointer to uint16_t positions: minus, plus, minus, plus, ...)
// pairs - a4  (number of pairs)
// out   - a5  (pointer to 12 int32_t sums, 16-byte aligned)
    entry   a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    ee.zero.q   q0          // Sums for shifts 0..3
    ee.zero.q   q1          // Sums for shifts 4..7
    ee.zero.q   q2          // Sums for shifts 8..11
    ee.zero.q   q7          // Zero for the sign comparison

    loopnez a4, .loop_end_sum_pairs_12_s16_pie
        l16ui   a6, a3, 0           // Position of the minus sample
        l16ui   a7, a3, 2           // Position of the plus sample
        addi    a3, a3, 4
        addx2   a6, a6, a2          // &in[minus]
        addx2   a7, a7, a2          // &in[plus]

        // Minus: samples 0..15 from &in[minus]
        ee.ld.128.usar.ip   q3, a6, 16
        ee.ld.128.usar.ip   q4, a6, 16
        ee.ld.128.usar.ip   q5, a6, 16
        ee.src.q    q3, q3, q4      // Samples 0..7
        ee.src.q    q4, q4, q5      // Samples 8..15
        ee.vcmp.lt.s16  q6, q3, q7  // High halves (sign) of samples 0..7
        ee.vzip.16  q3, q6          // q3 = samples 0..3, q6 = samples 4..7 (32 bit)
        ee.vsubs.s32    q0, q0, q3
        ee.vsubs.s32    q1, q1, q6
        ee.vcmp.lt.s16  q6, q4, q7  // High halves of samples 8..15
        ee.vzip.16  q4, q6          // q4 = samples 8..11 (32 bit)
        ee.vsubs.s32    q2, q2, q4

        // Plus: samples 0..15 from &in[plus]
        ee.ld.128.usar.ip   q3, a7, 16
        ee.ld.128.usar.ip   q4, a7, 16
        ee.ld.128.usar.ip   q5, a7, 16
        ee.src.q    q3, q3, q4
        ee.src.q    q4, q4, q5
        ee.vcmp.lt.s16  q6, q3, q7
        ee.vzip.16  q3, q6
        ee.vadds.s32    q0, q0, q3
        ee.vadds.s32    q1, q1, q6
        ee.vcmp.lt.s16  q6, q4, q7
        ee.vzip.16  q4, q6
        ee.vadds.s32    q2, q2, q4
    .loop_end_sum_pairs_12_s16_pie:

    ee.vst.128.ip   q0, a5, 16
    ee.vst.128.ip   q1, a5, 16
    ee.vst.128.ip   q2, a5, 16

    retw.n  // Return from function (windowed return)
