    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access)
    .text
    .type   magnitude_q15_pie,@function
    .align  4
    .global magnitude_q15_pie
// The function implements the following C code:
//void magnitude_q15(complex_q15* in, q15* out, uint32_t size)
//{
//    for(uint32_t i=0; i<size; i++)
//    {
//        int32_t x=in[i].re;
//        int32_t y=in[i].im;
//        out[i]=(((x*x) >> 16) + ((y*y) >> 16)); // Calculates magnitude squared divided by 2, in Q15 format
//    }
//}
magnitude_q15_pie: 
// in   - a2  (pointer to input complex array)
// out  - a3  (pointer to output magnitude array)
// size - a4  (number of complex numbers to process)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    srli 	a4, a4, 3  // Divide size by 8 (as the loop processes 8 elements per iteration)
    movi.n  a9, 1 	   // Load constant 1
    sub 	a4, a4, a9   // Decrement loop counter (a4 = size/8 - 1), used for loopnez
    ssai    16         // Set SAR register to 16 (used for vector multiplication shifts)
 
    ee.vld.128.ip           q0, a2, 16  // Load 128 bits (4 complex_q15 = 8 int16_t) from input into q0, increment a2 by 16
    ee.vmul.s16.ld.incp     q2, a2,q1, q0, q0   // q1 = q0 * q0 (square elements of q0), load next 128 bits into q2, increment a2 by 16 (q1 now holds Re^2, Im^2 for first 4 complex numbers)
    ee.vmul.s16             q3, q2, q2   // q3 = q2 * q2 (square elements of q2) (q3 now holds Re^2, Im^2 for next 4 complex numbers)
    ee.vunzip.16            q1, q3       // Interleave 16-bit elements of q1 and q3 (reorganize data for addition)
    ee.vadds.s16.ld.incp    q0, a2, q4, q1, q3   // q4 = q1 + q3 (add corresponding squared real and imaginary parts), load next 128 bits into q0, increment a2 by 16

    loopnez a4, .loop_end_magnitude_q15_pie  // Loop (size/8 - 1) times, exit if a4 becomes zero
        ee.vmul.s16.ld.incp    q2, a2,q1, q0, q0   // q1 = q0 * q0 (square elements), load next 128 bits into q2, increment a2 by 16
        ee.vmul.s16.st.incp    q4, a3, q3, q2, q2   // q3 = q2 * q2 (square elements), store previous result q4 to output, increment a3 by 16

        ee.vunzip.16           q1, q3               // Interleave 16-bit elements of q1 and q3
        ee.vadds.s16.ld.incp   q0, a2, q4, q1, q3   // q4 = q1 + q3 (add squared parts), load next 128 bits into q0, increment a2 by 16
.loop_end_magnitude_q15_pie:

    ee.vst.128.ip   q4, a3, 16  // Store the final calculated magnitude result (q4) to output array, increment a3 by 16
    retw.n  // Return from function (windowed return)
