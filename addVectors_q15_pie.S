
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
    entry	a1, 16

    srli 	a5, a5, 3 //  size/8 	

    ee.vld.128.ip   q0, a2, 16
    ee.vld.128.ip   q1, a3, 16
 
    loopnez a5, .loop_end_addVectors_q15_pie
        ee.vadds.s16.ld.incp    q0, a2, q2, q0, q1
        ee.vst.128.ip   q2, a4, 16
        ee.vld.128.ip   q1, a3, 16
.loop_end_addVectors_q15_pie:
    retw.n

