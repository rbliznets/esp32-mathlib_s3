    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access)
    .text
    .align  4
	.literal_position
	.literal .Q16_05, 2147450879  // Literal for 0.5 in Q16 format (used for normalization after addition/subtraction)
	.literal .Q15_05, 1073692671  // Literal for 0.5 in Q15 format (used for normalization after addition/subtraction)
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .text
    .type   fft_r2_q15_pie,@function
    .align  4
    .global fft_r2_q15_pie
fft_r2_q15_pie: 
// data  - a2  (pointer to input/output complex data array)
// w     - a3  (pointer to twiddle factor array)
// size  - a4  (FFT size)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    ssai    16 // Set SAR register to 16 (used for vector shifts later)

    l32r a7,.Q16_05     // Load 0.5 in Q16 format
    ee.movi.32.q	q7, a7, 0  // Load Q16_05 into all 4 32-bit slots of vector register q7
    ee.movi.32.q	q7, a7, 1 
    ee.movi.32.q	q7, a7, 2 
    ee.movi.32.q	q7, a7, 3
    l32r a6,.Q15_05     // Load 0.5 in Q15 format
    ee.movi.32.q	q6, a6, 0  // Load Q15_05 into all 4 32-bit slots of vector register q6
    ee.movi.32.q	q6, a6, 1 
    ee.movi.32.q	q6, a6, 2 
    ee.movi.32.q	q6, a6, 3
 
    movi.n a8,1     // s = 1 (stage size counter, number of butterflies per group in this stage)
    //movi.n a7,32    // (unused)
    slli a4,a4,2    // offset = size*4 (calculate byte offset for half the array, assuming complex_q15 is 4 bytes)
    
    add.n a6,a2,a4  // a6 = &data[size*4] (end address for outer loop comparison)
    addi.n a6,a6,16 // a6 = &data[size*4 + 16] (adjust end address)
    
.stageN:            // Stages with size > 4 (most general case)
    srai a4,a4,1        // offset = (offset/2) (halve the offset for the next stage)
    addi.n a12,a4,-16   // a12 = offset - 16 (used for incrementing the main data pointer in the loop)
 
    srai a13,a4,4       // a13 = (offset/16) (number of butterflies per group, adjusted for loop counter)
    mov.n a10,a2        // a10 = &data[0] (start of data for this stage)
 
.loop0_start_stageN:  // Outer loop: iterate through groups of butterflies for this stage
    add.n a11,a10,a4    // a11 = &data[i+offset] (pointer to second operand in butterfly)
    mov.n a14,a10       // a14 = &data[i] (save pointer for writing result back)
    mov.n a15,a11       // a15 = &data[i+offset] (save pointer for writing result back)
    mov.n a9,a3         // a9 = &w[0] (reset twiddle factor pointer for this group)

    EE.VLD.128.IP q1, a10, 16                   // q1 = data[i .. i+3] (load first 4 complex numbers), increment a10
    EE.VLD.128.IP q2, a11, 16                   // q2 = data[i+offset .. i+offset+3] (load second 4 complex numbers), increment a11
    EE.VADDS.S16.LD.INCP q0, a9, q3, q1, q2     // q0 = w[twiddle_index .. twiddle_index+3], q3 = data[i .. i+3] + data[i+offset .. i+offset+3]
    EE.VMUL.S16 q3, q3, q7                      // q3 = q3 * 0.5 (normalize the addition result)
    EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1    // data[i .. i+3] = q3 (store normalized sum), q4 = data[i+offset .. i+offset+3] - data[i .. i+3]
    EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[i+4 .. i+7], q4 = q4 * w[twiddle_index .. twiddle_index+3] (multiply difference by twiddle)
    EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[i+offset+4 .. i+offset+7], q4 = q4 * w[twiddle_index+4 .. twiddle_index+7] (multiply difference by next twiddle)

    addi.n a5,a13,-1    // Calculate loop count for inner butterflies (a13 - 1)
    loopnez a5, .loop1_end_stageN  // Inner loop: process remaining butterflies in the group
        EE.VADDS.S16.LD.INCP q0, a9, q3, q1, q2  // q0 = w[next_twiddle_index .. next_twiddle_index+3], q3 = data[i+4 .. i+7] + data[i+offset+4 .. i+offset+7]
        EE.VMUL.S16.ST.INCP q4, a15, q3, q3, q7  // data[i+offset+4 .. i+offset+7] = q4 (store previous difference result), q3 = q3 * 0.5 (normalize new sum)
        EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1 // data[i+4 .. i+7] = q3 (store normalized sum), q4 = data[i+offset+4 .. i+offset+7] - data[i+4 .. i+7] (calculate new difference)
        EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[next 4], q4 = q4 * w (multiply difference by twiddle)
        EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[next 4], q4 = q4 * w (multiply difference by next twiddle)
    .loop1_end_stageN:
    EE.VST.128.IP q4, a15, 16                   // data[i+offset+last_4] = q4 (store final difference result for this group), increment a15
    add.n a10,a10,a12   // Move to the next group of butterflies (a10 += offset - 16)
    bne a6, a11, .loop0_start_stageN  // Continue outer loop if we haven't processed all groups
    slli a8,a8,1        // s *= 2 (double the stage size for the next iteration)
    mov.n a3,a9         // Update twiddle factor pointer (a9 was advanced during the loop)
    bnei a4,32,.stageN  // Continue stages N if offset is not 32 (i.e., not the 3rd stage from the end)

.stage3:                // Stage where butterfly size is 4 (special case)
    mov.n a10,a2        // a10 = &data[0] (start of data for this stage)
    addi.n a11,a10,16   // a11 = &data[offset] (pointer to second operand, offset=4*4=16 bytes for 4 complex numbers)
    mov.n a14,a10       // a14 = &data[i] (save pointer for writing result back)
    mov.n a15,a11       // a15 = &data[i+offset] (save pointer for writing result back)

    EE.VLD.128.IP q1, a10, 32                   // q1 = data[i .. i+3] (load first 4 complex numbers), increment a10 by 32 (skip next 4 for now)
    EE.VLD.128.IP q2, a11, 32                   // q2 = data[i+offset .. i+offset+3] (load second 4 complex numbers), increment a11 by 32
    EE.VADDS.S16.LD.INCP q0, a9, q3, q1, q2     // q0 = w[twiddle_index .. twiddle_index+3], q3 = data[i .. i+3] + data[i+offset .. i+offset+3]
    EE.VMUL.S16 q3, q3, q7                      // q3 = q3 * 0.5 (normalize the addition result)
    EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1    // data[i .. i+3] = q3 (store normalized sum), q4 = data[i+offset .. i+offset+3] - data[i .. i+3]
    EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[i+8 .. i+11] (reload from new address), q4 = q4 * w (multiply difference by twiddle)
    EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[i+offset+8 .. i+offset+11] (reload from new address), q4 = q4 * w (multiply difference by next twiddle)
 
    addi.n a5,a8,-1     // Calculate number of groups for this stage (a8 butterflies per group, so a8-1 iterations)
    loopnez a5, .loop1_end_stage3  // Inner loop: process remaining groups for this stage
        addi.n a10,a10,16 // Move data pointers by 16 bytes (4 complex numbers)
        addi.n a11,a11,16

        EE.VADDS.S16.LD.INCP q5, a14, q3, q1, q2    // a14+=16, q3 = data[next 4] + data[next 4] (from updated a14)
        EE.VMUL.S16.ST.INCP q4, a15, q3, q3, q7     // data[next 4 offset] = q4 (store previous difference), q3 = q3 * 0.5 (normalize new sum)
        EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1    // data[next 4] = q3 (store normalized sum), q4 = data[next 4 offset] - data[next 4] (calculate new difference)
        EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[next 4], q4 = q4 * w (multiply difference by twiddle)
        EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[next 4], q4 = q4 * w (multiply difference by next twiddle)
        
        addi.n a15,a15,16 // Move output pointer by 16 bytes
    .loop1_end_stage3:
    EE.VST.128.IP q4, a15, 16                   // data[final_offset] = q4 (store final difference result), increment a15
 
 .stage2:               // Stage where butterfly size is 2 (special case)
    ssai    15      // Set SAR register to 15 (used for shifts in this stage)
    mov.n a10,a2    // a10 = &data[0] (start of data for this stage)
    mov.n a11,a2    // a11 = &data[0] (start of data for output in this stage)

    EE.VLD.H.64.IP q0, a9, 0                    // q0hi = w[0 1] (load first two twiddle factors into high part of q0)
    EE.VLD.128.IP q1, a10, 16                   // q1 = data[i i+1 i+offset i+offset+1] (load 2 complex pairs)
    EE.VLD.128.IP q2, a10, 16                   // q2 = data[i+2 i+3 i+offset+2 i+offset+3] (load next 2 complex pairs)
    EE.FFT.R2BF.S16 q5, q7, q1, q2, 0           // Perform radix-2 butterfly operation: q5 = out1, q7 = out2
    EE.VMUL.S16 q3, q5, q6                      // q3 = q5 * 0.5 (normalize out1)
    EE.VMUL.S16 q4, q7, q6                      // q4 = q7 * 0.5 (normalize out2)
    EE.CMUL.S16.LD.INCP q1, a10, q3, q3, q0, 1  // q3hi = q3hi * w[0 1] (multiply high parts of normalized result by twiddle factors)
    EE.CMUL.S16.LD.INCP q2, a10, q4, q4, q0, 1  // q4hi = q4hi * w[0 1] (multiply high parts of normalized result by twiddle factors)
    EE.VST.128.IP q3, a11, 16                   // data[i i+1 i+offset i+offset+1] = q3 (store first result pair), increment a11
    addi.n a5,a8,-1     // Calculate number of groups for this stage (a8 butterflies per group, so a8-1 iterations)
    loopnez a5, .loop1_end_stage2  // Inner loop: process remaining groups for this stage
        EE.FFT.R2BF.S16 q5, q7, q1, q2, 0           // Perform radix-2 butterfly operation
        EE.VMUL.S16.ST.INCP q4, a11, q3, q5, q6    // data[next pair offset] = q4 (store previous result), q3 = q5 * 0.5 (normalize new out1)
        EE.VMUL.S16.LD.INCP q1, a10, q4, q7, q6    // q1 = data[next pair], q4 = q7 * 0.5 (normalize new out2)

        EE.CMUL.S16.LD.INCP q2, a10, q3, q3, q0, 1 // q2 = data[next pair], q3hi = q3hi * w (multiply normalized out1 by twiddle)
        EE.CMUL.S16.ST.INCP q3, a11, q4, q4, q0, 1 // data[next pair] = q3 (store result from previous iteration), q4hi = q4hi * w (multiply normalized out2 by twiddle)
    .loop1_end_stage2:
    EE.VST.128.IP q4, a11, 16 // data[final pair] = q4 (store final result), increment a11

.stage1:                // Stage where butterfly size is 1 (final stage, trivial butterflies with w=1)
    mov.n a10,a2        // a10 = &data[0] (start of data for this stage)
    mov.n a11,a2        // a11 = &data[0] (start of data for output in this stage)
    EE.VLDBC.32 q0, a9  // q0 = w[0 0 0 0] (broadcast first twiddle factor (should be 1+j0) to all slots)
    EE.VLD.128.IP q1, a10, 16 // q1 = data[i i+1 i+offset i+offset+1] (load first pair of complex numbers)
    EE.VLD.128.IP q2, a10, 16 // q2 = data[i+2 i+3 i+offset+2 i+offset+3] (load next pair of complex numbers)
    loopnez a8, .loop1_end_stage1  // Loop a8 times (number of groups in the final stage)
        EE.FFT.R2BF.S16 q3, q4, q1, q2, 1  // Perform radix-2 butterfly operation (final stage logic, likely just add/sub)

        EE.CMUL.S16.LD.INCP q1, a10, q3, q3, q0, 1  // q1 = data[next pair], q3hi = q3hi * w[0] (multiply result by final twiddle)
        EE.CMUL.S16.LD.INCP q2, a10, q4, q4, q0, 1  // q2 = data[next pair], q4hi = q4hi * w[0] (multiply result by final twiddle)
        EE.VUNZIP.32 q3, q4               // Interleave 32-bit elements of q3 and q4 (reorganize results)
        
        EE.VST.L.64.IP q3, a11, 8  // Store low 64 bits of q3, increment a11 by 8
        EE.VST.L.64.IP q4, a11, 8  // Store low 64 bits of q4, increment a11 by 8
        EE.VST.H.64.IP q3, a11, 8  // Store high 64 bits of q3, increment a11 by 8
        EE.VST.H.64.IP q4, a11, 8  // Store high 64 bits of q4, increment a11 by 8
    .loop1_end_stage1:

    retw.n  // Return from function (windowed return)
