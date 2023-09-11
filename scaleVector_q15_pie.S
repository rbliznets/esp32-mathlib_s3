// The function implements the following C code:
//void scaleVector(q15* in, q15* k, q15* out, uint32_t size)
//{
//    for(uint32_t i=0; i<size; i++)
//    {
//        out[i]=(in1[i]*(*k)) >> 15;
//    }
//}

    .text
    .type   scaleVector_q15_pie,@function
    .align  4
    .global scaleVector_q15_pie
scaleVector_q15_pie: 
// in   - a2
// k   - a3
// out   - a4
// size  - a5
    entry	a1, 16

    srli 	a5, a5, 3 //  size/8 
    movi.n  a6, 1
    sub     a5, a5, a6	// size--

    EE.LD.128.USAR.IP   q3, a2, 16 
    ee.vldbc.16     q1, a3     // q1 = 8[k]
    ssai            15         // sar = 15 

    EE.LD.128.USAR.IP   q4, a2, 16 
    EE.SRC.Q.QUP        q0, q3, q4
    EE.VMUL.S16         q2, q0, q1

    loopnez a5, .loop_end_scaleVector_q15_pie
        EE.LD.128.USAR.IP   q4, a2, 16 
        EE.SRC.Q.QUP        q0, q3, q4
        EE.VMUL.S16.ST.INCP q2, a4, q2, q0, q1
    .loop_end_scaleVector_q15_pie:
    ee.vst.128.ip       q2, a4, 16
    retw.n
/////////////////////////////////////////////////////////////////////////////////////////////////
    .type   scaleVector_q15_16_pie,@function
    .align  4
    .global scaleVector_q15_16_pie
scaleVector_q15_16_pie: 
// in   - a2
// k   - a3
// out   - a4
// size  - a5
    entry	a1, 16

    srli 	a5, a5, 3 //  size/8 	

    ee.vld.128.ip   q0, a2, 16
    ee.vldbc.16     q1, a3     // q1 = 8[k]
    ssai            15         // sar = 15 

    loopnez a5, .loop_end_scaleVector_q15_16_pie
        ee.vmul.s16.ld.incp     q0, a2, q2, q0, q1
        ee.vst.128.ip           q2, a4, 16
    .loop_end_scaleVector_q15_16_pie:
    retw.n
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
// in   - a2
// shift   - a3
// out   - a4
// size  - a5
    entry	a1, 16

    srli 	a5, a5, 3 //  size/8 
    WSR.SAR a3

    loopnez a5, .loop_end_shrinkVector_16_pie
        ee.vld.128.ip   q0,a2,16
        ee.vld.128.ip   q1,a2,16
        ee.vsr.32       q0,q0
        ee.vsr.32       q1,q1
        ee.vunzip.16    q0,q1
        ee.vst.128.ip   q0, a4, 16
    .loop_end_shrinkVector_16_pie:
    retw.n

