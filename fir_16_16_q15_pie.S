    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access)
    .text
    .align  4
	.literal_position
	.literal	.SAT_MAX_30, 1073692671   // Literal for maximum value that fits in 30 bits after accumulation and saturation (around 2^29.999 - 1)
	.literal	.SAT_MIN_30, -1073692672  // Literal for minimum value that fits in 30 bits after accumulation and saturation (around -2^29.999)
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   fir_16_16_q15_pie,@function
    .align  4
    .global fir_16_16_q15_pie
fir_16_16_q15_pie: 
// in    - a2  (pointer to input data array, 16-byte aligned)
// k     - a3  (pointer to FIR coefficients array, 16-byte aligned)
// ksize - a4  (number of coefficients in the FIR filter, must be multiple of 8)
// out   - a5  (pointer to output data array, 16-byte aligned)
// size  - a6  (number of output samples to generate, must be multiple of 8)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    srli 	a4, a4, 3           // Divide ksize by 8 (as the inner loop processes 8 coefficients per iteration)
    slli 	a6, a6, 1           // Multiply size by 2 (calculate byte offset for end of output array)
    add.n 	a6, a6, a5          // a6 = out + size*2 (calculate end address of output array)

    EE.VLD.128.IP  q0, a2, 16   // q0 = in[0 .. 7] (load first 8 input samples), increment a2 by 16
    MOV.N a7, a2                // a7 = &in[8] (save pointer to next input block, advanced by 8 samples)
    MOV.N a8, a3                // a8 = k (save the initial coefficient pointer)
    MOVI.N  a9, 15              // a9 = 15 (used later as the shift amount for normalization from internal 30-bit to Q15)

.big_loop:                      // Outer loop: process one output sample per iteration (though the code handles 8 samples per loop iteration due to vectorization)
    EE.ZERO.QACC                        // Clear the Q accumulator (QACC = 0) for accumulating the dot product
    loopnez a4, .loop_end_fir_pie       // Inner loop: iterate ksize/8 times
        EE.VLD.128.IP  q1, a2, 16       // q1 = in[i+8 .. i+15] (load next 8 input samples), increment a2 by 16
        EE.VLD.128.IP  q5, a3, 16       // q5 = k[i .. i+7] (load 8 coefficients), increment a3 by 16

        EE.VSMULAS.S16.QACC q0, q5, 0   // QACC[0..7] += in[i .. i+7] * k[i] (multiply-add 8 pairs, accumulate in QACC)
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+1 .. i+8], q1 = in[i+9 .. i+15 0] (shift q0 left by 1 element, shifting in from q1; q1 shifts left, filling low with 0)
        EE.VSMULAS.S16.QACC q0, q5, 1   // QACC[0..7] += in[i+1 .. i+8] * k[i+1] (multiply-add next 8 pairs)
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+2 .. i+9], q1 = in[i+10 .. i+15 0 0] (shift again)
        EE.VSMULAS.S16.QACC q0, q5, 2   // QACC[0..7] += in[i+2 .. i+9] * k[i+2]
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+3 .. i+10], q1 = in[i+11 .. i+15, 0 .. 0] (shift again)
        EE.VSMULAS.S16.QACC q0, q5, 3   // QACC[0..7] += in[i+3 .. i+10] * k[i+3]
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+4 .. i+11], q1 = in[i+12 .. i+15, 0 .. 0] (shift again)
        EE.VSMULAS.S16.QACC q0, q5, 4   // QACC[0..7] += in[i+4 .. i+11] * k[i+4]
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+5 .. i+12], q1 = in[i+13 .. i+15, 0 .. 0] (shift again)
        EE.VSMULAS.S16.QACC q0, q5, 5   // QACC[0..7] += in[i+5 .. i+12] * k[i+5]
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+6 .. i+13], q1 = in[i+14 i+15, 0 .. 0] (shift again)
        EE.VSMULAS.S16.QACC q0, q5, 6   // QACC[0..7] += in[i+6 .. i+13] * k[i+6]
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+7 .. i+14], q1 = in[i+15, 0 .. 0] (shift again)
        EE.VSMULAS.S16.QACC q0, q5, 7   // QACC[0..7] += in[i+7 .. i+14] * k[i+7]
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+8 .. i+15], q1 = 0 (shift again, q1 becomes all zeros)
    .loop_end_fir_pie:                  // End of inner loop (processed one set of 8 coefficients)
    EE.SRCMB.S16.QACC q1, a9, 0         // q1 = saturate(QACC[0..7] >> 15) (saturate accumulated result and shift right by 15 to normalize from internal 30-bit to Q15 format)
    EE.VST.128.IP q1, a5, 16            // out[i..i+7] = q1 (store 8 output samples), increment a5 by 16
    BEQ a5,a6,.end_fir_16_16_q15_pie    // Check if output pointer a5 has reached the end a6; if so, jump to end

    EE.VLD.128.IP  q0, a7, 16       // q0 = in[i .. i+7] (load next block of 8 input samples from saved pointer a7), increment a7 by 16
    MOV.N a2, a7                    // a2 = a7 (update main input pointer to the next block)
    MOV.N a3, a8                    // a3 = a8 (reset coefficient pointer to the beginning)
    J .big_loop                     // Jump back to the start of the outer loop

.end_fir_16_16_q15_pie:
    retw.n  // Return from function (windowed return)
