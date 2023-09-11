
    .text
    .align  4
	.literal_position
	.literal .LC1_MX_32, 2147483647
	.literal .Q15_05, 1073692671
/////////////////////////////////////////////////////////////////////////////////////////
    .type   normalize_q15_pie,@function
    .align  4
    .global normalize_q15_pie
normalize_q15_pie: 
// in   - a2
// out   - a3
// size  - a4
    entry	a1, 16

    srli 	a4, a4, 3 //  size/8
    mov     a5, a2 	
    mov     a8, a4 	

    ee.vld.128.ip   q0, a2, 16
    ee.zero.q       q1 //minimum
    ee.zero.q       q2 //maximum
 
    loopnez a4, .loop1_end_normalize_q15_pie
        ee.vmin.s16             q1, q0, q1          //minimum x8
        ee.vmax.s16.ld.incp     q0, a2, q2, q0, q2  //maximum x8
.loop1_end_normalize_q15_pie:
    ee.vunzip.16    q1, q3
    ee.vmin.s16     q1, q3, q1 //minimum x4
    ee.vunzip.16    q1, q3
    ee.vmin.s16     q1, q3, q1 //minimum x2
    ee.vunzip.16    q1, q3
    ee.vmin.s16     q1, q3, q1 //minimum x1

    ee.vunzip.16    q2, q3
    ee.vmax.s16     q2, q3, q2 //maximum x4
    ee.vunzip.16    q2, q3
    ee.vmax.s16     q2, q3, q2 //maximum x2
    ee.vunzip.16    q2, q3
    ee.vmax.s16     q2, q3, q2 //maximum x1

    ee.movi.32.a	q1, a6, 0 //minimum
    slli 	        a6, a6, 16 
    abs 	        a6, a6 
    ee.movi.32.a	q2, a7, 0 //maximum
    slli 	        a7, a7, 16 
    max             a2, a6, a7 //max(-minimum,maximum)

    bnez    a2, .l1_normalize_q15_pie 
        retw.n  //return if max == 0
.l1_normalize_q15_pie:

    l32r	a6, .LC1_MX_32 
    quos    a2, a6, a2  // mult=0x7fffffff/(max << 16)
    bnei    a2, 1, .l2_normalize_q15_pie
        retw.n //return if mult == 1
.l2_normalize_q15_pie:

    slli            a6, a2, 16
    or              a6, a2, a6
    ee.movi.32.q	q1, a6, 0 
    ee.movi.32.q	q1, a6, 1 
    ee.movi.32.q	q1, a6, 2 
    ee.movi.32.q	q1, a6, 3
    ssai            0         // sar = 0 
    ee.vld.128.ip   q0, a5, 16
    loopnez a8, .loop2_end_normalize_q15_pie
        ee.vmul.s16.ld.incp     q0, a5, q2, q0, q1  //scale
        ee.vst.128.ip   q2, a3, 16
.loop2_end_normalize_q15_pie:

    retw.n
/////////////////////////////////////////////////////////////////////////////////////////
    .type   normalize_q14_pie,@function
    .align  4
    .global normalize_q14_pie
normalize_q14_pie: 
// in   - a2
// out   - a3
// size  - a4
    entry	a1, 16

    srli 	a4, a4, 3 //  size/8
    mov     a5, a2 	
    mov     a8, a4 	

    ee.vld.128.ip   q0, a2, 16
    ee.zero.q       q1 //minimum
    ee.zero.q       q2 //maximum
 
    loopnez a4, .loop1_end_normalize_q14_pie
        ee.vmin.s16             q1, q0, q1          //minimum x8
        ee.vmax.s16.ld.incp     q0, a2, q2, q0, q2  //maximum x8
.loop1_end_normalize_q14_pie:
    ee.vunzip.16    q1, q3
    ee.vmin.s16     q1, q3, q1 //minimum x4
    ee.vunzip.16    q1, q3
    ee.vmin.s16     q1, q3, q1 //minimum x2
    ee.vunzip.16    q1, q3
    ee.vmin.s16     q1, q3, q1 //minimum x1

    ee.vunzip.16    q2, q3
    ee.vmax.s16     q2, q3, q2 //maximum x4
    ee.vunzip.16    q2, q3
    ee.vmax.s16     q2, q3, q2 //maximum x2
    ee.vunzip.16    q2, q3
    ee.vmax.s16     q2, q3, q2 //maximum x1

    ee.movi.32.a	q1, a6, 0 //minimum
    slli 	        a6, a6, 16 
    abs 	        a6, a6 
    ee.movi.32.a	q2, a7, 0 //maximum
    slli 	        a7, a7, 16 
    max             a2, a6, a7 //max(-minimum,maximum)

    bnez    a2, .l1_normalize_q14_pie 
        retw.n  //return if max == 0
.l1_normalize_q14_pie:

    l32r	a6, .LC1_MX_32 
    quos    a2, a6, a2  // mult=0x7fffffff/(max << 16)
    bnei    a2, 1, .l2_normalize_q14_pie
    l32r a14,.Q15_05
    ee.movi.32.q	q1, a14, 0 
    ee.movi.32.q	q1, a14, 1 
    ee.movi.32.q	q1, a14, 2 
    ee.movi.32.q	q1, a14, 3
    ssai            15         // sar = 15 
    j .normalize_loop
.l2_normalize_q14_pie:
    blti a2,4,.loop2_end_normalize_q14_pie
    slli            a6, a2, 16
    or              a6, a2, a6
    ee.movi.32.q	q1, a6, 0 
    ee.movi.32.q	q1, a6, 1 
    ee.movi.32.q	q1, a6, 2 
    ee.movi.32.q	q1, a6, 3
    ssai            1         // sar = 1
 .normalize_loop:
    ee.vld.128.ip   q0, a5, 16
    loopnez a8, .loop2_end_normalize_q14_pie
        ee.vmul.s16.ld.incp     q0, a5, q2, q0, q1  //scale
        ee.vst.128.ip   q2, a3, 16
.loop2_end_normalize_q14_pie:

    retw.n
