    .text
    .align  4
	.literal_position
	.literal .Q16_05, 2147450879
	.literal .Q15_05, 1073692671
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .text
    .type   fft_r2_q15_pie,@function
    .align  4
    .global fft_r2_q15_pie
fft_r2_q15_pie: 
// data  - a2
// w     - a3
// size  - a4
    entry	a1, 16

    ssai    16 // sar = 16

    l32r a7,.Q16_05
    ee.movi.32.q	q7, a7, 0 
    ee.movi.32.q	q7, a7, 1 
    ee.movi.32.q	q7, a7, 2 
    ee.movi.32.q	q7, a7, 3
    l32r a6,.Q15_05
    ee.movi.32.q	q6, a6, 0 
    ee.movi.32.q	q6, a6, 1 
    ee.movi.32.q	q6, a6, 2 
    ee.movi.32.q	q6, a6, 3
 
    movi.n a8,1     // s = 1
    //movi.n a7,32    // 
    slli a4,a4,2    // offset = size*4 (bytes)
    
    add.n a6,a2,a4
    addi.n a6,a6,16
.stageN:
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
    EE.VLD.128.IP q2, a11, 16                   // q2 = data[i+offset .. i+offset+3]
    EE.VADDS.S16.LD.INCP q0, a9, q3, q1, q2     //  q0 = w[i .. i+3], q3 =data[i .. i+3]+data[i+offset .. i+offset+3]
    EE.VMUL.S16 q3, q3, q7                      // q3 = norm*q3
    EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1    // data[i .. i+3] = q3, q4 = data[i+offset .. i+offset+3] - data[i .. i+3]
    EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[i .. i+3], q4 = q4 * w
    EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[i+offset .. i+offset+3], q4 = q4 * w

    addi.n a5,a13,-1
    loopnez a5, .loop1_end_stageN
        EE.VADDS.S16.LD.INCP q0, a9, q3, q1, q2  // q0 = w[i .. i+3], q3 =data[i .. i+3]+data[i+offset .. i+offset+3]
        EE.VMUL.S16.ST.INCP q4, a15, q3, q3, q7  // data[i+offset .. i+offset+3] = q4, q3 = norm*q3
        EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1 // data[i .. i+3] = q3, q4 = data[i+offset .. i+offset+3] - data[i .. i+3]
        EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[i .. i+3], q4 = q4 * w
        EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[i+offset .. i+offset+3], q4 = q4 * w
    .loop1_end_stageN:
    EE.VST.128.IP q4, a15, 16                   //  data[i+offset .. i+offset+3] = q4
    add.n a10,a10,a12
    bne a6, a11, .loop0_start_stageN
    slli a8,a8,1
    mov.n a3,a9
    bnei a4,32,.stageN

.stage3:
    mov.n a10,a2        // a10 = &data[0] (load)
    addi.n a11,a10,16   // a11 = &data[offset] (load)
    mov.n a14,a10       // a14 = &data[i] (save)
    mov.n a15,a11       // a15 = &data[i+offset] (save)

    EE.VLD.128.IP q1, a10, 32                   // q1 = data[i .. i+3]
    EE.VLD.128.IP q2, a11, 32                   // q2 = data[i+offset .. i+offset+3]
    EE.VADDS.S16.LD.INCP q0, a9, q3, q1, q2     //  q0 = w[i .. i+3], q3 =data[i .. i+3]+data[i+offset .. i+offset+3]
    EE.VMUL.S16 q3, q3, q7                      // q3 = norm*q3
    EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1    // data[i .. i+3] = q3, q4 = data[i+offset .. i+offset+3] - data[i .. i+3]
    EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[i .. i+3], q4 = q4 * w
    EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[i+offset .. i+offset+3], q4 = q4 * w
 
    addi.n a5,a8,-1
    loopnez a5, .loop1_end_stage3
        addi.n a10,a10,16
        addi.n a11,a11,16

        EE.VADDS.S16.LD.INCP q5, a14, q3, q1, q2    // a14+=16, q3 =data[i .. i+3]+data[i+offset .. i+offset+3]
        EE.VMUL.S16.ST.INCP q4, a15, q3, q3, q7     // data[i+offset .. i+offset+3] = q4, q3 = norm*q3
        EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1    // data[i .. i+3] = q3, q4 = data[i+offset .. i+offset+3] - data[i .. i+3]
        EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[i .. i+3], q4 = q4 * w
        EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[i+offset .. i+offset+3], q4 = q4 * w
        
        addi.n a15,a15,16
    .loop1_end_stage3:
    EE.VST.128.IP q4, a15, 16                   //  data[i+offset .. i+offset+3] = q4
 
 .stage2:
    ssai    15      // sar = 15
    mov.n a10,a2    // a10 = &data[0] (load)
    mov.n a11,a2    // a11 = &data[0] (save)

    EE.VLD.H.64.IP q0, a9, 0                    // q0hi = w[0 1]
    EE.VLD.128.IP q1, a10, 16                   // q1 = data[i i+1 i+offset i+offset+1]
    EE.VLD.128.IP q2, a10, 16                   // q2 = data[i+2 i+3 i+offset+2 i+offset+3]
    EE.FFT.R2BF.S16 q5, q7, q1, q2, 0           // butterfly
    EE.VMUL.S16 q3, q5, q6                      // q3 = norm*q5
    EE.VMUL.S16 q4, q7, q6                      // q4 = norm*q7
    EE.CMUL.S16.LD.INCP q1, a10, q3, q3, q0, 1  //q3hi = q3hi * w
    EE.CMUL.S16.LD.INCP q2, a10, q4, q4, q0, 1  // q4hi = q4hi * w
    EE.VST.128.IP q3, a11, 16                   // data[i i+1 i+offset i+offset+1] = q3
    addi.n a5,a8,-1
    loopnez a5, .loop1_end_stage2
        EE.FFT.R2BF.S16 q5, q7, q1, q2, 0           // butterfly
        EE.VMUL.S16.ST.INCP q4, a11, q3, q5, q6    // data[i+2 i+3 i+offset+2 i+offset+3] = q4 ,q3 = norm*q5
        EE.VMUL.S16.LD.INCP q1, a10, q4, q7, q6    // q1 = data[i i+1 i+offset i+offset+1], q4 = norm*q7

        EE.CMUL.S16.LD.INCP q2, a10, q3, q3, q0, 1 // q1 = data[i i+1 i+offset i+offset+1], q3hi = q3hi * w
        EE.CMUL.S16.ST.INCP q3, a11, q4, q4, q0, 1 // data[i i+1 i+offset i+offset+1] = q3, q4hi = q4hi * w
    .loop1_end_stage2:
    EE.VST.128.IP q4, a11, 16 // data[i+2 i+3 i+offset+2 i+offset+3] = q4

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

    retw.n
