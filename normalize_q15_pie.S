    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access)
    .text
    .align  4
	.literal_position
	.literal .LC1_MX_32, 2147483647   // Literal for maximum value of a 32-bit signed integer (0x7FFFFFFF)
	.literal .Q15_05, 1073692671      // Literal for 0.5 in Q15 format (used when mult=1)
/////////////////////////////////////////////////////////////////////////////////////////
    .type   normalize_q15_pie,@function
    .align  4
    .global normalize_q15_pie
normalize_q15_pie: 
// in   - a2  (pointer to input array)
// out  - a3  (pointer to output array)
// size - a4  (number of elements to process, multiple of 8)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    srli 	a4, a4, 3 // Divide size by 8 (as the loop processes 8 elements per iteration)
    mov     a5, a2 	   // Save the initial input pointer to a5
    mov     a8, a4 	   // Save the loop count (size/8) to a8
 
    ee.vld.128.ip   q0, a2, 16  // Load 128 bits (8 int16_t) from input into q0, increment a2 by 16
    ee.zero.q       q1 // Initialize minimum value accumulator q1 to zero
    ee.zero.q       q2 // Initialize maximum value accumulator q2 to zero
 
    loopnez a4, .loop1_end_normalize_q15_pie  // Loop size/8 times
        ee.vmin.s16             q1, q0, q1          // Update minimum values in q1 using elements from q0
        ee.vmax.s16.ld.incp     q0, a2, q2, q0, q2  // Update maximum values in q2 using elements from q0, load next 128 bits into q0, increment a2 by 16
.loop1_end_normalize_q15_pie:
    // Find the single minimum value across all elements in q1
    ee.vunzip.16    q1, q3
    ee.vmin.s16     q1, q3, q1 //minimum x4 (reduce within pairs)
    ee.vunzip.16    q1, q3
    ee.vmin.s16     q1, q3, q1 //minimum x2 (reduce further)
    ee.vunzip.16    q1, q3
    ee.vmin.s16     q1, q3, q1 //minimum x1 (reduce to single value in q1)

    // Find the single maximum value across all elements in q2
    ee.vunzip.16    q2, q3
    ee.vmax.s16     q2, q3, q2 //maximum x4 (reduce within pairs)
    ee.vunzip.16    q2, q3
    ee.vmax.s16     q2, q3, q2 //maximum x2 (reduce further)
    ee.vunzip.16    q2, q3
    ee.vmax.s16     q2, q3, q2 //maximum x1 (reduce to single value in q2)

    ee.movi.32.a	q1, a6, 0 // Move the minimum value (low 32 bits of q1) to register a6
    slli 	        a6, a6, 16  // Shift minimum left by 16 bits (to get the high 16 bits of the int16_t)
    abs 	        a6, a6      // Take absolute value of the minimum
    ee.movi.32.a	q2, a7, 0 // Move the maximum value (low 32 bits of q2) to register a7
    slli 	        a7, a7, 16  // Shift maximum left by 16 bits
    max             a2, a6, a7 // a2 = max(abs(minimum), maximum) (find the largest absolute value in the data)

    bnez    a2, .l1_normalize_q15_pie  // If the maximum absolute value is not zero, continue normalization
        retw.n  // Return if max == 0 (no normalization needed, all values are zero)
.l1_normalize_q15_pie:

    l32r	a6, .LC1_MX_32  // Load the maximum 32-bit signed integer value (0x7FFFFFFF)
    quos    a2, a6, a2  // a2 = 0x7FFFFFFF / (max_abs_value << 16) (calculate scaling factor)
    bnei    a2, 1, .l2_normalize_q15_pie  // If scaling factor is not 1, continue with general normalization
        retw.n // Return if mult == 1 (no scaling needed, values already fit in Q15 range)
.l2_normalize_q15_pie:

    slli            a6, a2, 16    // a6 = scaling_factor << 16
    or              a6, a2, a6    // a6 = scaling_factor | (scaling_factor << 16) (create 16-bit repeating pattern)
    ee.movi.32.q	q1, a6, 0     // Load the scaling factor into all slots of vector register q1
    ee.movi.32.q	q1, a6, 1 
    ee.movi.32.q	q1, a6, 2 
    ee.movi.32.q	q1, a6, 3
    ssai            0         // Set SAR register to 0 (for the normalization multiplication shift)
    ee.vld.128.ip   q0, a5, 16  // Load first 128 bits of input data into q0, increment a5 by 16
    loopnez a8, .loop2_end_normalize_q15_pie  // Loop size/8 times
        ee.vmul.s16.ld.incp     q0, a5, q2, q0, q1  // Multiply q0 by scaling factor q1, store result in q2, load next data into q0, increment a5 by 16 (scale)
        ee.vst.128.ip   q2, a3, 16  // Store the scaled result from q2 to output array, increment a3 by 16
.loop2_end_normalize_q15_pie:

    retw.n  // Return from function (windowed return)
/////////////////////////////////////////////////////////////////////////////////////////
    .type   normalize_q14_pie,@function
    .align  4
    .global normalize_q14_pie
normalize_q14_pie: 
// in   - a2  (pointer to input array)
// out  - a3  (pointer to output array)
// size - a4  (number of elements to process, multiple of 8)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    srli 	a4, a4, 3 // Divide size by 8 (as the loop processes 8 elements per iteration)
    mov     a5, a2 	   // Save the initial input pointer to a5
    mov     a8, a4 	   // Save the loop count (size/8) to a8
 
    ee.vld.128.ip   q0, a2, 16  // Load 128 bits (8 int16_t) from input into q0, increment a2 by 16
    ee.zero.q       q1 // Initialize minimum value accumulator q1 to zero
    ee.zero.q       q2 // Initialize maximum value accumulator q2 to zero
 
    loopnez a4, .loop1_end_normalize_q14_pie  // Loop size/8 times
        ee.vmin.s16             q1, q0, q1          // Update minimum values in q1 using elements from q0
        ee.vmax.s16.ld.incp     q0, a2, q2, q0, q2  // Update maximum values in q2 using elements from q0, load next 128 bits into q0, increment a2 by 16
.loop1_end_normalize_q14_pie:
    // Find the single minimum value across all elements in q1 (similar to q15 version)
    ee.vunzip.16    q1, q3
    ee.vmin.s16     q1, q3, q1 //minimum x4
    ee.vunzip.16    q1, q3
    ee.vmin.s16     q1, q3, q1 //minimum x2
    ee.vunzip.16    q1, q3
    ee.vmin.s16     q1, q3, q1 //minimum x1

    // Find the single maximum value across all elements in q2 (similar to q15 version)
    ee.vunzip.16    q2, q3
    ee.vmax.s16     q2, q3, q2 //maximum x4
    ee.vunzip.16    q2, q3
    ee.vmax.s16     q2, q3, q2 //maximum x2
    ee.vunzip.16    q2, q3
    ee.vmax.s16     q2, q3, q2 //maximum x1

    ee.movi.32.a	q1, a6, 0 // Move the minimum value (low 32 bits of q1) to register a6
    slli 	        a6, a6, 16 
    abs 	        a6, a6 
    ee.movi.32.a	q2, a7, 0 // Move the maximum value (low 32 bits of q2) to register a7
    slli 	        a7, a7, 16 
    max             a2, a6, a7 // a2 = max(abs(minimum), maximum) (find the largest absolute value in the data)

    bnez    a2, .l1_normalize_q14_pie  // If the maximum absolute value is not zero, continue normalization
        retw.n  // Return if max == 0 (no normalization needed, all values are zero)
.l1_normalize_q14_pie:

    l32r	a6, .LC1_MX_32  // Load the maximum 32-bit signed integer value (0x7FFFFFFF)
    quos    a2, a6, a2  // a2 = 0x7FFFFFFF / (max_abs_value << 16) (calculate scaling factor)
    bnei    a2, 1, .l2_normalize_q14_pie  // If scaling factor is not 1, continue with general normalization
    // Scaling factor is 1, but for Q14, we still need to scale by 0.5 (represented as Q15_05) to fit in the Q14 range
    l32r a14,.Q15_05      // Load 0.5 in Q15 format
    ee.movi.32.q	q1, a14, 0  // Load 0.5 into all slots of q1
    ee.movi.32.q	q1, a14, 1 
    ee.movi.32.q	q1, a14, 2 
    ee.movi.32.q	q1, a14, 3
    ssai            15         // Set SAR register to 15 (for Q15 multiplication shift)
    j .normalize_loop        // Jump to the normalization loop
.l2_normalize_q14_pie:
    blti a2,4,.loop2_end_normalize_q14_pie  // If scaling factor is less than 4, skip normalization (values are small enough for Q14)
    slli            a6, a2, 16    // a6 = scaling_factor << 16
    or              a6, a2, a6    // a6 = scaling_factor | (scaling_factor << 16) (create 16-bit repeating pattern)
    ee.movi.32.q	q1, a6, 0     // Load the scaling factor into all slots of vector register q1
    ee.movi.32.q	q1, a6, 1 
    ee.movi.32.q	q1, a6, 2 
    ee.movi.32.q	q1, a6, 3
    ssai            1         // Set SAR register to 1 (for the normalization multiplication shift)
 .normalize_loop:
    ee.vld.128.ip   q0, a5, 16  // Load first 128 bits of input data into q0, increment a5 by 16
    loopnez a8, .loop2_end_normalize_q14_pie  // Loop size/8 times
        ee.vmul.s16.ld.incp     q0, a5, q2, q0, q1  // Multiply q0 by scaling factor q1, store result in q2, load next data into q0, increment a5 by 16 (scale)
        ee.vst.128.ip   q2, a3, 16  // Store the scaled result from q2 to output array, increment a3 by 16
.loop2_end_normalize_q14_pie:

    retw.n  // Return from function (windowed return)

