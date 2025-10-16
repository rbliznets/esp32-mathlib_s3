    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access)
    .text
    .align  4
	.literal_position
	.literal	.SAT_MAX_30, 1073741823    // Literal for maximum value that fits in 30 bits after accumulation (2^29 - 1)
	.literal	.SAT_MIN_30, -1073741824   // Literal for minimum value that fits in 30 bits after accumulation (-2^29)
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   dot_product_q15_16_16_pie,@function
    .align  4
    .global dot_product_q15_16_16_pie
dot_product_q15_16_16_pie: 
// in1   - a2  (pointer to first input vector, 16-byte aligned)
// in2   - a3  (pointer to second input vector, 16-byte aligned)
// size  - a4  (number of elements to multiply and sum)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    srli 	a4, a4, 3           // Divide size by 8 (as the loop processes 8 elements per iteration)
    ee.zero.accx                // Clear the accumulator register (accx = 0)
    ee.vld.128.ip   q0, a2, 16  // Load 128 bits (8 int16_t) from in1[0..7] into q0, increment a2 by 16
    ee.vld.128.ip   q1, a3, 16  // Load 128 bits (8 int16_t) from in2[0..7] into q1, increment a3 by 16

    loopnez a4, .loop_end_dot_product_q15_16_16_pie  // Loop 'size/8' times, exit if a4 becomes zero
        ee.vmulas.s16.accx.ld.ip    q0, a2, 16, q0, q1  // Multiply-add signed 16-bit elements (q0 * q1), accumulate in accx, load next 128 bits from in1 into q0, increment a2 by 16
        ee.vld.128.ip               q1, a3, 16          // Load next 128 bits (8 int16_t) from in2 into q1, increment a3 by 16
    .loop_end_dot_product_q15_16_16_pie:

    movi.n	a5, 0
    EE.SRS.ACCX a2, a5, 0  // Saturate the accumulator (accx) and store the result in a2

    l32r	a3, .SAT_MAX_30    // Load the maximum 30-bit saturated value
    min     a2, a2, a3          // Clamp result to maximum if overflow occurred during saturation
    l32r	a4, .SAT_MIN_30    // Load the minimum 30-bit saturated value
    max     a2, a2, a4          // Clamp result to minimum if underflow occurred during saturation
    SRAI    a2, a2, 15          // Arithmetic right shift by 15 (divide by 2^15 to normalize from internal 30-bit format to Q15)

    retw.n  // Return from function (windowed return)
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   dot_product_q15_1_16_pie,@function
    .align  4
    .global dot_product_q15_1_16_pie
dot_product_q15_1_16_pie: 
// in1   - a2  (pointer to first input vector, potentially unaligned)
// in2   - a3  (pointer to second input vector, 16-byte aligned)
// size  - a4  (number of elements to multiply and sum)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    srli 	a4, a4, 3           // Divide size by 8 (as the loop processes 8 elements per iteration)

    ee.zero.accx                // Clear the accumulator register (accx = 0)
    EE.LD.128.USAR.IP   q2, a2, 16  // Load 128 bits from in1 using SAR shift, increment a2 by 16 (handles unaligned access)
    ee.vld.128.ip   q1, a3, 16  // Load 128 bits (8 int16_t) from in2[0..7] into q1, increment a3 by 16

    loopnez a4, .loop_end_dot_product_q15_1_16_pie  // Loop 'size/8' times, exit if a4 becomes zero
        EE.LD.128.USAR.IP   q3, a2, 16  // Load next 128 bits from in1 using SAR shift, increment a2 by 16
        EE.SRC.Q.QUP        q0, q2, q3  // Shift and combine q2 and q3 to reorganize data into q0 (e.g., unpacking or aligning elements from the unaligned load)
        ee.vmulas.s16.accx.ld.ip    q1, a3, 16, q0, q1  // Multiply-add signed 16-bit elements (q0 * q1), accumulate in accx, load next 128 bits from in2 into q1, increment a3 by 16
    .loop_end_dot_product_q15_1_16_pie:

    movi.n	a5, 0
    EE.SRS.ACCX a2, a5, 0  // Saturate the accumulator (accx) and store the result in a2

    l32r	a3, .SAT_MAX_30    // Load the maximum 30-bit saturated value
    min     a2, a2, a3          // Clamp result to maximum if overflow occurred during saturation
    l32r	a4, .SAT_MIN_30    // Load the minimum 30-bit saturated value
    max     a2, a2, a4          // Clamp result to minimum if underflow occurred during saturation
    SRAI    a2, a2, 15          // Arithmetic right shift by 15 (divide by 2^15 to normalize from internal 30-bit format to Q15)

    retw.n  // Return from function (windowed return)
