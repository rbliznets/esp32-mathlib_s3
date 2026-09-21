    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access)
    .align  4
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Bit-reversal permutation of complex_q15 data in place for FFT sizes 16..1024.
//
// EE.BITREV qa, as writes into qa eight 16-bit swap partners for the indices as..as+7 and adds 8 to as.
// The partner of x is bitrev(x) (FFT_BIT_WIDTH bits) if it is greater than x, otherwise x itself, so
// every index is swapped with its partner unconditionally: each pair is exchanged once (by its smaller
// index) and a swap with itself changes nothing. FFT_BIT_WIDTH holds the number of bits modulo 8
// (8..10 bits -> 0..2); 3 bits (size 8) give wrong partners, so sizes below 16 are not supported.
// Two lanes are swapped together: both loads go before both stores, which is safe because the partner
// of x is never x+1 (bitrev(x) = x+1 would need an even x below N/2 with the top bit set).
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   revbin_permute_pie,@function
    .align  4
    .global revbin_permute_pie
revbin_permute_pie:
// data  - a2  (pointer to complex_q15 array, 4-byte aligned)
// size  - a3  (FFT size, power of 2, 16..1024)
// width - a4  (log2(size) & 7 for the FFT_BIT_WIDTH register)
    entry   a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    wur.fft_bit_width a4    // Number of reversed bits (modulo 8)
    ssai    2               // SAR = 2: EE.VSL.32 turns an index into a byte offset
    ee.movi.32.q    q6, a2, 0
    ee.movi.32.q    q6, a2, 1
    ee.movi.32.q    q6, a2, 2
    ee.movi.32.q    q6, a2, 3   // q6 = data address in all four 32-bit lanes
    movi.n  a5, 0           // Index counter for EE.BITREV (incremented by 8 by the instruction)
    srli    a6, a3, 3       // Number of groups of 8 indices
    mov.n   a7, a2          // a7 = &data[x] for the first index of the group

    loopnez a6, .loop_end_revbin_permute_pie
        ee.bitrev       q0, a5          // q0 = partners of x..x+7 (16 bit), a5 += 8
        ee.zero.q       q1
        ee.vzip.16      q0, q1          // zero-extend: q0 = partners 0..3, q1 = partners 4..7 (32 bit)
        ee.vsl.32       q0, q0          // byte offsets
        ee.vsl.32       q1, q1
        ee.vadds.s32    q0, q0, q6      // addresses of the partners
        ee.vadds.s32    q1, q1, q6

        ee.movi.32.a    q0, a9, 0       // &data[partner] of lane
        ee.movi.32.a    q0, a8, 1
        l32i    a10, a7, 0
        l32i    a11, a9, 0
        l32i    a12, a7, 4
        l32i    a13, a8, 0
        s32i    a11, a7, 0
        s32i    a10, a9, 0
        s32i    a13, a7, 4
        s32i    a12, a8, 0

        ee.movi.32.a    q0, a9, 2       // &data[partner] of lane
        ee.movi.32.a    q0, a8, 3
        l32i    a10, a7, 8
        l32i    a11, a9, 0
        l32i    a12, a7, 12
        l32i    a13, a8, 0
        s32i    a11, a7, 8
        s32i    a10, a9, 0
        s32i    a13, a7, 12
        s32i    a12, a8, 0

        ee.movi.32.a    q1, a9, 0       // &data[partner] of lane
        ee.movi.32.a    q1, a8, 1
        l32i    a10, a7, 16
        l32i    a11, a9, 0
        l32i    a12, a7, 20
        l32i    a13, a8, 0
        s32i    a11, a7, 16
        s32i    a10, a9, 0
        s32i    a13, a7, 20
        s32i    a12, a8, 0

        ee.movi.32.a    q1, a9, 2       // &data[partner] of lane
        ee.movi.32.a    q1, a8, 3
        l32i    a10, a7, 24
        l32i    a11, a9, 0
        l32i    a12, a7, 28
        l32i    a13, a8, 0
        s32i    a11, a7, 24
        s32i    a10, a9, 0
        s32i    a13, a7, 28
        s32i    a12, a8, 0

        addi    a7, a7, 32              // next group
    .loop_end_revbin_permute_pie:

    retw.n  // Return from function (windowed return)
