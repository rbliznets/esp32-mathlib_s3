    .section .iram1,"ax"
    .align  4
	.literal_position
	.literal .Q15_05, 1073692671
	.literal .LC1_MX_32, 2147483647
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   fft_r2s_q15_pie,@function
    .align  4
    .global fft_r2s_q15_pie
fft_r2s_q15_pie: 
// data  - a2
// w     - a3
// size  - a4
    entry	a1, 16

    ssai    15 // sar = 15

    l32r a6,.Q15_05
    ee.movi.32.q	q5, a6, 0 
    ee.movi.32.q	q5, a6, 1 
    ee.movi.32.q	q5, a6, 2 
    ee.movi.32.q	q5, a6, 3
 
    movi.n a8,1     // s = 1
    srai a7,a4,2    // 
    slli a4,a4,2    // offset = size*4 (bytes)
    
    add.n a6,a2,a4
    addi.n a6,a6,16
.stageN:
    ee.zero.q q6 //minimum
    ee.zero.q q7 //maximum

    srai a4,a4,1        // offset = (offset/2)
    addi.n a12,a4,-16
 
    srai a13,a4,4       // (offset/8)
    mov.n a10,a2        // a10 = &data[0]
 
.loop0_start_stageN:
    add.n a11,a10,a4    // a11 = &data[i+offset] (load)
    mov.n a14,a10       // a14 = &data[i] (save)
    mov.n a15,a11       // a15 = &data[i+offset] (save)
    mov.n a9,a3         // a9 = &w[0]

    EE.VLD.128.IP q1, a10, 16                   // q1 = data[i .. i+3]
    EE.VLD.128.IP q2, a11, 16      // q2 = data[i+offset .. i+offset+3], q1 = norm*q1
    EE.VADDS.S16.LD.INCP q0, a9, q3, q1, q2     //  q0 = w[i .. i+3], q3 =data[i .. i+3]+data[i+offset .. i+offset+3]
    EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1    // data[i .. i+3] = q3, q4 = data[i+offset .. i+offset+3] - data[i .. i+3]
    EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[i .. i+3], q4 = q4 * w
    EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[i+offset .. i+offset+3], q4 = q4 * w
    ee.vmin.s16.ST.INCP q4, a15, q6, q3, q6  //minimum x8
    ee.vmax.s16  q7, q3, q7  //maximum x8
    ee.vmin.s16  q6, q4, q6  //minimum x8
    ee.vmax.s16  q7, q4, q7  //maximum x8

    addi.n a5,a13,-1
    loopnez a5, .loop1_end_stageN
        EE.VADDS.S16.LD.INCP q0, a9, q3, q1, q2  // q0 = w[i .. i+3], q3 =data[i .. i+3]+data[i+offset .. i+offset+3]
        EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1 // data[i .. i+3] = q3, q4 = data[i+offset .. i+offset+3] - data[i .. i+3]
        EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[i .. i+3], q4 = q4 * w
        EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[i+offset .. i+offset+3], q4 = q4 * w
        ee.vmin.s16.ST.INCP q4, a15,  q6, q3, q6  //minimum x8
        ee.vmax.s16  q7, q3, q7  //maximum x8
        ee.vmin.s16  q6, q4, q6  //minimum x8
        ee.vmax.s16  q7, q4, q7  //maximum x8
   .loop1_end_stageN:
    add.n a10,a10,a12
    bne a6, a11, .loop0_start_stageN

    ee.vunzip.16    q6, q4
    ee.vmin.s16     q6, q4, q6 //minimum x4
    ee.vunzip.16    q6, q4
    ee.vmin.s16     q6, q4, q6 //minimum x2
    ee.vunzip.16    q6, q4
    ee.vmin.s16     q6, q4, q6 //minimum x1
    ee.vunzip.16    q7, q4
    ee.vmax.s16     q7, q4, q7 //maximum x4
    ee.vunzip.16    q7, q4
    ee.vmax.s16     q7, q4, q7 //maximum x2
    ee.vunzip.16    q7, q4
    ee.vmax.s16     q7, q4, q7 //maximum x1
    ee.movi.32.a	q6, a5, 0 //minimum
    slli 	        a5, a5, 16 
    abs 	        a5, a5
    ee.movi.32.a	q7, a14, 0 //maximum
    slli 	        a14, a14, 16 
    max             a5, a5, a14 //max(-minimum,maximum)

    beqz    a5, .end 
    l32r	a15, .LC1_MX_32 
    quos    a5, a15, a5  // mult=0x7fffffff/(max << 16)
    bnei    a5,1,.normalize_stageN_0 
    l32r a14,.Q15_05
    ee.movi.32.q	q5, a14, 0 
    ee.movi.32.q	q5, a14, 1 
    ee.movi.32.q	q5, a14, 2 
    ee.movi.32.q	q5, a14, 3
    j .normalize_stageN
.normalize_stageN_0: 
    blti a5,4,.end_normalize_stageN
    slli            a15, a5, 16
    or              a15, a5, a15
    ee.movi.32.q	q5, a15, 0 
    ee.movi.32.q	q5, a15, 1 
    ee.movi.32.q	q5, a15, 2 
    ee.movi.32.q	q5, a15, 3
    ssai            1         // sar = 1 
.normalize_stageN:
    mov.n a10,a2    // a10 = &data[0] (load)
    mov.n a11,a2    // a11 = &data[0] (save)
    ee.vld.128.ip   q1, a10, 16
    loopnez a7, .stageN_end_normalize
        ee.vmul.s16.ld.incp     q1, a10, q2, q1, q5  //scale
        ee.vst.128.ip   q2, a11, 16
    .stageN_end_normalize:
    ssai            15         // sar = 15 
.end_normalize_stageN:

    slli a8,a8,1
    mov.n a3,a9
    bnei a4,32,.stageN

.stage3:
    ee.zero.q q6 //minimum
    ee.zero.q q7 //maximum

    mov.n a10,a2        // a10 = &data[0] (load)
    addi.n a11,a10,16   // a11 = &data[offset] (load)
    mov.n a14,a10       // a14 = &data[i] (save)
    mov.n a15,a11       // a15 = &data[i+offset] (save)

    EE.VLD.128.IP q1, a10, 32                   // q1 = data[i .. i+3]
    EE.VLD.128.IP q2, a11, 32                   // q2 = data[i+offset .. i+offset+3]
    EE.VADDS.S16.LD.INCP q0, a9, q3, q1, q2     //  q0 = w[i .. i+3], q3 =data[i .. i+3]+data[i+offset .. i+offset+3]
    EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1    // data[i .. i+3] = q3, q4 = data[i+offset .. i+offset+3] - data[i .. i+3]
    EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[i .. i+3], q4 = q4 * w
    EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[i+offset .. i+offset+3], q4 = q4 * w
    ee.vmin.s16.ST.INCP q4, a15,  q6, q3, q6  //minimum x8
    ee.vmax.s16.LD.INCP q5, a10,    q7, q3, q7  //maximum x8
    ee.vmin.s16.LD.INCP q5, a11,    q6, q4, q6  //minimum x8
    ee.vmax.s16.LD.INCP q5, a14,    q7, q4, q7  //maximum x8
 
    addi.n a5,a8,-1
    loopnez a5, .loop1_end_stage3
        EE.VADDS.S16.LD.INCP q5, a15, q3, q1, q2    // a14+=16, q3 =data[i .. i+3]+data[i+offset .. i+offset+3]
        EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1    // data[i .. i+3] = q3, q4 = data[i+offset .. i+offset+3] - data[i .. i+3]
        EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[i .. i+3], q4 = q4 * w
        EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[i+offset .. i+offset+3], q4 = q4 * w
        
        ee.vmin.s16.ST.INCP q4, a15,  q6, q3, q6  //minimum x8
        ee.vmax.s16.LD.INCP q5, a10,    q7, q3, q7  //maximum x8
        ee.vmin.s16.LD.INCP q5, a11,    q6, q4, q6  //minimum x8
        ee.vmax.s16.LD.INCP q5, a14,    q7, q4, q7  //maximum x8
    .loop1_end_stage3:

    ee.vunzip.16    q6, q4
    ee.vmin.s16     q6, q4, q6 //minimum x4
    ee.vunzip.16    q6, q4
    ee.vmin.s16     q6, q4, q6 //minimum x2
    ee.vunzip.16    q6, q4
    ee.vmin.s16     q6, q4, q6 //minimum x1
    ee.vunzip.16    q7, q4
    ee.vmax.s16     q7, q4, q7 //maximum x4
    ee.vunzip.16    q7, q4
    ee.vmax.s16     q7, q4, q7 //maximum x2
    ee.vunzip.16    q7, q4
    ee.vmax.s16     q7, q4, q7 //maximum x1
    ee.movi.32.a	q6, a5, 0 //minimum
    slli 	        a5, a5, 16 
    abs 	        a5, a5
    ee.movi.32.a	q7, a14, 0 //maximum
    slli 	        a14, a14, 16 
    max             a5, a5, a14 //max(-minimum,maximum)

    beqz    a5, .end 
    l32r	a15, .LC1_MX_32 
    quos    a5, a15, a5  // mult=0x7fffffff/(max << 16)
    bnei    a5,1,.normalize_stage2_0 
    l32r a14,.Q15_05
    ee.movi.32.q	q5, a14, 0 
    ee.movi.32.q	q5, a14, 1 
    ee.movi.32.q	q5, a14, 2 
    ee.movi.32.q	q5, a14, 3
    j .normalize_stage2
.normalize_stage2_0: 
    blti a5,4,.stage2
    slli            a15, a5, 16
    or              a15, a5, a15
    ee.movi.32.q	q5, a15, 0 
    ee.movi.32.q	q5, a15, 1 
    ee.movi.32.q	q5, a15, 2 
    ee.movi.32.q	q5, a15, 3
    ssai            1         // sar = 1 
.normalize_stage2:
    mov.n a10,a2    // a10 = &data[0] (load)
    mov.n a11,a2    // a11 = &data[0] (save)
    ee.vld.128.ip   q1, a10, 16
    loopnez a7, .stage2_end_normalize
        ee.vmul.s16.ld.incp     q1, a10, q2, q1, q5  //scale
        ee.vst.128.ip   q2, a11, 16
    .stage2_end_normalize:
    ssai            15         // sar = 15 
////////////////////////////////////////////////////////////////////////// 
 .stage2:
    ee.zero.q q6 //minimum
    ee.zero.q q7 //maximum

    mov.n a10,a2    // a10 = &data[0] (load)
    mov.n a11,a2    // a11 = &data[0] (save)

    EE.VLD.H.64.IP q0, a9, 0                    // q0hi = w[0 1]
    EE.VLD.128.IP q1, a10, 16                   // q1 = data[i i+1 i+offset i+offset+1]
    EE.VLD.128.IP q2, a10, 16     //  q2 = data[i+2 i+3 i+offset+2 i+offset+3]
    loopnez a8, .loop1_end_stage2
        EE.FFT.R2BF.S16 q3, q4, q1, q2, 0  // butterfly

        EE.CMUL.S16.LD.INCP q1, a10, q3, q3, q0, 1  //q3hi = q3hi * w
        EE.CMUL.S16.LD.INCP q2, a10, q4, q4, q0, 1 // q4hi = q4hi * w

        ee.vmin.s16.ST.INCP q3, a11, q6, q3, q6  //minimum x8
        ee.vmax.s16.ST.INCP q4, a11, q7, q3, q7  //maximum x8
        ee.vmin.s16  q6, q4, q6  //minimum x8
        ee.vmax.s16  q7, q4, q7  //maximum x8
    .loop1_end_stage2:
////////////////////////////////////
    ee.vunzip.16    q6, q4
    ee.vmin.s16     q6, q4, q6 //minimum x4
    ee.vunzip.16    q6, q4
    ee.vmin.s16     q6, q4, q6 //minimum x2
    ee.vunzip.16    q6, q4
    ee.vmin.s16     q6, q4, q6 //minimum x1
    ee.vunzip.16    q7, q4
    ee.vmax.s16     q7, q4, q7 //maximum x4
    ee.vunzip.16    q7, q4
    ee.vmax.s16     q7, q4, q7 //maximum x2
    ee.vunzip.16    q7, q4
    ee.vmax.s16     q7, q4, q7 //maximum x1
    ee.movi.32.a	q6, a5, 0 //minimum
    slli 	        a5, a5, 16 
    abs 	        a5, a5
    ee.movi.32.a	q7, a14, 0 //maximum
    slli 	        a14, a14, 16 
    max             a5, a5, a14 //max(-minimum,maximum)

//@@@@@@@@@@@@@@@@@@@@@@@@@@
    //S32I.N a5, a2, 0
    //retw.n


    beqz    a5, .end 
    l32r	a15, .LC1_MX_32 
    quos    a5, a15, a5  // mult=0x7fffffff/(max << 16)
    bnei    a5,1,.normalize_stage1_0 
    l32r a14,.Q15_05
    ee.movi.32.q	q5, a14, 0 
    ee.movi.32.q	q5, a14, 1 
    ee.movi.32.q	q5, a14, 2 
    ee.movi.32.q	q5, a14, 3
    j .normalize_stage1
.normalize_stage1_0: 
    blti a5,4,.stage1
    slli            a15, a5, 16
    or              a15, a5, a15
    ee.movi.32.q	q5, a15, 0 
    ee.movi.32.q	q5, a15, 1 
    ee.movi.32.q	q5, a15, 2 
    ee.movi.32.q	q5, a15, 3
    ssai            1         // sar = 1 
.normalize_stage1:
    mov.n a10,a2    // a10 = &data[0] (load)
    mov.n a11,a2    // a11 = &data[0] (save)
    ee.vld.128.ip   q1, a10, 16
    loopnez a7, .stage1_end_normalize
        ee.vmul.s16.ld.incp     q1, a10, q2, q1, q5  //scale
        ee.vst.128.ip   q2, a11, 16
    .stage1_end_normalize:
    ssai            15         // sar = 15 
////////////////////////////////////
.stage1:
    mov.n a10,a2        // a10 = &data[0] (load)
    mov.n a11,a2        // a11 = &data[0] (save)
    EE.VLDBC.32 q0, a9  // q0 = w[0 0 0 0]
    EE.VLD.128.IP q1, a10, 16 // q1 = data[i i+1 i+offset i+offset+1]
    EE.VLD.128.IP q2, a10, 16 // q2 = data[i+2 i+3 i+offset+2 i+offset+3]
    loopnez a8, .loop1_end_stage1
        EE.FFT.R2BF.S16 q3, q4, q1, q2, 1  // butterfly

        EE.CMUL.S16.LD.INCP q1, a10, q3, q3, q0, 1  //q3hi = q3hi * w
        EE.CMUL.S16.LD.INCP q2, a10, q4, q4, q0, 1  // q4hi = q4hi * w
        EE.VUNZIP.32 q3, q4
        
        EE.VST.L.64.IP q3, a11, 8 
        EE.VST.L.64.IP q4, a11, 8 
        EE.VST.H.64.IP q3, a11, 8 
        EE.VST.H.64.IP q4, a11, 8 
    .loop1_end_stage1:

.end:
    retw.n
