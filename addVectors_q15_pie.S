    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access)
    .text
    .type   addVectors_q15_pie,@function
    .align  4
    .global addVectors_q15_pie
// The function implements the following C code:
//void addVectors(int16_t* in1, int16_t* in2, int16_t* out, uint32_t size)
//{
//    for(uint32_t i=0; i<size; i++)
//    {
//        out[i]=in1[i]+in2[i];
//    }
//}
addVectors_q15_pie: 
// in1   - a2
// in2   - a3
// out   - a4
// size  - a5
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    srli 	a5, a5, 3 // Divide size by 8 (since we process 8 elements per loop iteration)

    ee.vld.128.ip   q0, a2, 16  // Load 128 bits (8 x int16_t) from in1 into vector register q0, increment in1 pointer by 16 bytes
    ee.vld.128.ip   q1, a3, 16  // Load 128 bits (8 x int16_t) from in2 into vector register q1, increment in2 pointer by 16 bytes
 
    loopnez a5, .loop_end_addVectors_q15_pie  // Loop 'size/8' times, exit if a5 becomes zero
        ee.vadds.s16.ld.incp    q0, a2, q2, q0, q1  // Add signed 16-bit integers from q0 and q1, saturate result into q2, load next 128 bits from updated in1 pointer into q0, increment in1 pointer again
        ee.vst.128.ip   q2, a4, 16               // Store 128 bits (8 x int16_t) from q2 to out, increment out pointer by 16 bytes
        ee.vld.128.ip   q1, a3, 16               // Load 128 bits (8 x int16_t) from updated in2 pointer into q1, increment in2 pointer by 16 bytes
.loop_end_addVectors_q15_pie:
    retw.n  // Return from function (windowed return)
