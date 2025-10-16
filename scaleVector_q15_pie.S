// The function implements the following C code:
//void scaleVector(q15* in, q15* k, q15* out, uint32_t size)
//{
//    for(uint32_t i=0; i<size; i++)
//    {
//        out[i]=(in[i]*(*k)) >> 15; // Multiply each input element by a scalar and normalize by shifting right by 15
//    }
//}

    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access)
    .text
    .type   scaleVector_q15_pie,@function
    .align  4
    .global scaleVector_q15_pie
scaleVector_q15_pie: 
// in   - a2  (pointer to input array, potentially unaligned)
// k    - a3  (pointer to scalar value)
// out  - a4  (pointer to output array, 16-byte aligned)
// size - a5  (number of elements to process, must be >= 16)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    srli 	a5, a5, 3 // Divide size by 8 (as the loop processes 8 elements per iteration)
    movi.n  a6, 1
    sub     a5, a5, a6	// Decrement loop counter (a5 = size/8 - 1), used for loopnez

    EE.LD.128.USAR.IP   q3, a2, 16  // Load 128 bits from input (a2) using SAR shift, increment a2 by 16 (handles unaligned access)
    ee.vldbc.16     q1, a3     // q1 = broadcast 8 copies of the scalar value pointed to by a3
    ssai            15         // Set SAR register to 15 (for the multiplication shift)

    EE.LD.128.USAR.IP   q4, a2, 16  // Load next 128 bits from input using SAR shift, increment a2 by 16
    EE.SRC.Q.QUP        q0, q3, q4  // Shift and combine q3 and q4 to reorganize data into q0 (aligns elements from unaligned loads)
    EE.VMUL.S16         q2, q0, q1  // q2 = q0 * q1 (multiply aligned input elements by the scalar)

    loopnez a5, .loop_end_scaleVector_q15_pie  // Loop (size/8 - 1) times
        EE.LD.128.USAR.IP   q4, a2, 16  // Load next 128 bits from input using SAR shift, increment a2 by 16
        EE.SRC.Q.QUP        q0, q3, q4  // Shift and combine to reorganize data into q0
        EE.VMUL.S16.ST.INCP q2, a4, q2, q0, q1  // q2 = q0 * q1 (multiply), store previous result q2 to output, increment a4 by 16
    .loop_end_scaleVector_q15_pie:
    ee.vst.128.ip       q2, a4, 16  // Store the final calculated result q2 to output array, increment a4 by 16
    retw.n  // Return from function (windowed return)
/////////////////////////////////////////////////////////////////////////////////////////////////
    .type   scaleVector_q15_16_pie,@function
    .align  4
    .global scaleVector_q15_16_pie
scaleVector_q15_16_pie: 
// in   - a2  (pointer to input array, 16-byte aligned)
// k    - a3  (pointer to scalar value, 2-byte aligned)
// out  - a4  (pointer to output array, 16-byte aligned)
// size - a5  (number of elements to process, multiple of 8)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    srli 	a5, a5, 3 // Divide size by 8 (as the loop processes 8 elements per iteration)

    ee.vld.128.ip   q0, a2, 16  // Load 128 bits (8 int16_t) from input into q0, increment a2 by 16
    ee.vldbc.16     q1, a3     // q1 = broadcast 8 copies of the scalar value pointed to by a3
    ssai            15         // Set SAR register to 15 (for the multiplication shift)

    loopnez a5, .loop_end_scaleVector_q15_16_pie  // Loop size/8 times
        ee.vmul.s16.ld.incp     q0, a2, q2, q0, q1  // q2 = q0 * q1 (multiply), load next 128 bits into q0, increment a2 by 16
        ee.vst.128.ip           q2, a4, 16          // Store the result q2 to output array, increment a4 by 16
    .loop_end_scaleVector_q15_16_pie:
    retw.n  // Return from function (windowed return)
/////////////////////////////////////////////////////////////////////////////////////////////////
//	complex_q15* cdata = (complex_q15*)data;
//	for (int i = 0; i < (size/2); i++)
//	{
//		cdata[i].re = (q15)(data[2*i] >> mLastShift);
//		cdata[i].im = (q15)(data[2*i+1] >> mLastShift);
//	}
    .type   shrinkVector_16_pie,@function
    .align  4
    .global shrinkVector_16_pie
shrinkVector_16_pie: 
// in    - a2  (pointer to input array of uint32_t, 16-byte aligned)
// shift - a3  (number of bits to shift right, 1-31)
// out   - a4  (pointer to output array of q15, 16-byte aligned)
// size  - a5  (number of elements to process, multiple of 8)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    srli 	a5, a5, 3 // Divide size by 8 (as the loop processes 8 elements per iteration)
    WSR.SAR a3        // Write the shift amount (a3) to the SAR (Shift Amount Register)

    loopnez a5, .loop_end_shrinkVector_16_pie  // Loop size/8 times
        ee.vld.128.ip   q0,a2,16  // Load 128 bits (4 uint32_t) from input into q0, increment a2 by 16
        ee.vld.128.ip   q1,a2,16  // Load next 128 bits (4 uint32_t) from input into q1, increment a2 by 16
        ee.vsr.32       q0,q0     // Shift each 32-bit element in q0 right by the SAR amount (arithmetic shift)
        ee.vsr.32       q1,q1     // Shift each 32-bit element in q1 right by the SAR amount (arithmetic shift)
        ee.vunzip.16    q0,q1     // Interleave 16-bit elements of q0 and q1 (e.g., q0 = [q0_low, q1_low], q1 = [q0_high, q1_high])
        ee.vst.128.ip   q0, a4, 16  // Store the reorganized and shifted 16-bit values from q0 to output, increment a4 by 16
    .loop_end_shrinkVector_16_pie:
    retw.n  // Return from function (windowed return)
