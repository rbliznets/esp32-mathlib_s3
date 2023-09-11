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
//        out[i]=(((x*x) >> 16) + ((y*y) >> 16));
//    }
//}
magnitude_q15_pie: 
// in   - a2
// out   - a3
// size   - a4
    entry	a1, 16

    srli 	a4, a4, 3  //  size/8
    movi.n  a9, 1 	
    sub 	a4, a4, a9  	
    ssai    16         // sar = 16
 
    ee.vld.128.ip           q0, a2, 16
    ee.vmul.s16.ld.incp     q2, a2,q1, q0, q0   //q1 = a^2,b^2
    ee.vmul.s16             q3, q2, q2   //q3 = a^2,b^2
    ee.vunzip.16            q1, q3
    ee.vadds.s16.ld.incp    q0, a2, q4, q1, q3   //q4=a^2+b^2

    loopnez a4, .loop_end_magnitude_q15_pie
        ee.vmul.s16.ld.incp    q2, a2,q1, q0, q0   //q1 = a^2,b^2
        ee.vmul.s16.st.incp    q4, a3, q3, q2, q2   //q3 = a^2,b^2

        ee.vunzip.16           q1, q3
        ee.vadds.s16.ld.incp   q0, a2, q4, q1, q3   //q4=a^2+b^2
.loop_end_magnitude_q15_pie:

    ee.vst.128.ip   q4, a3, 16
    retw.n

