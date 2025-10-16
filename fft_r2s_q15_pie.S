    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access)
    .align  4
	.literal_position
	.literal .Q15_05, 1073692671    // Literal for 0.5 in Q15 format (used for normalization)
	.literal .LC1_MX_32, 2147483647  // Literal for maximum value of a 32-bit signed integer (0x7FFFFFFF)
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   fft_r2s_q15_pie,@function
    .align  4
    .global fft_r2s_q15_pie
fft_r2s_q15_pie: 
// data  - a2  (pointer to input/output complex data array)
// w     - a3  (pointer to twiddle factor array)
// size  - a4  (FFT size)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    ssai    15 // Set SAR register to 15 (used for vector shifts later)

    l32r a6,.Q15_05     // Load 0.5 in Q15 format
    ee.movi.32.q	q5, a6, 0  // Load Q15_05 into all 4 32-bit slots of vector register q5 (used for scaling)
    ee.movi.32.q	q5, a6, 1 
    ee.movi.32.q	q5, a6, 2 
    ee.movi.32.q	q5, a6, 3
 
    movi.n a8,1     // s = 1 (stage size counter, number of butterflies per group in this stage)
    srai a7,a4,2    // a7 = size / 4 (used later for normalization loop count)
    slli a4,a4,2    // offset = size*4 (calculate byte offset for half the array, assuming complex_q15 is 4 bytes)
    
    add.n a6,a2,a4  // a6 = &data[size*4] (end address for outer loop comparison)
    addi.n a6,a6,16 // a6 = &data[size*4 + 16] (adjust end address)
.stageN:            // Stages with size > 4 (most general case)
    ee.zero.q q6 // Initialize minimum value accumulator q6 to zero
    ee.zero.q q7 // Initialize maximum value accumulator q7 to zero

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
    EE.VLD.128.IP q2, a11, 16      // q2 = data[i+offset .. i+offset+3] (load second 4 complex numbers), q1 = norm*q1 (likely a typo in original comment)
    EE.VADDS.S16.LD.INCP q0, a9, q3, q1, q2     // q0 = w[i .. i+3], q3 = data[i .. i+3] + data[i+offset .. i+offset+3]
    EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1    // data[i .. i+3] = q3 (store sum), q4 = data[i+offset .. i+offset+3] - data[i .. i+3] (calculate difference)
    EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[i+4 .. i+7], q4 = q4 * w (multiply difference by twiddle)
    EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[i+offset+4 .. i+offset+7], q4 = q4 * w (multiply difference by next twiddle)
    ee.vmin.s16.ST.INCP q4, a15, q6, q3, q6  // data[i+offset+4 .. i+offset+7] = q4 (store difference), update minimum in q6 using q3 values
    ee.vmax.s16  q7, q3, q7  // Update maximum in q7 using q3 values
    ee.vmin.s16  q6, q4, q6  // Update minimum in q6 using q4 values
    ee.vmax.s16  q7, q4, q7  // Update maximum in q7 using q4 values

    addi.n a5,a13,-1    // Calculate loop count for inner butterflies (a13 - 1)
    loopnez a5, .loop1_end_stageN  // Inner loop: process remaining butterflies in the group
        EE.VADDS.S16.LD.INCP q0, a9, q3, q1, q2  // q0 = w[next_twiddle_index .. next_twiddle_index+3], q3 = data[i+4 .. i+7] + data[i+offset+4 .. i+offset+7]
        EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1 // data[i+4 .. i+7] = q3 (store sum), q4 = data[i+offset+4 .. i+offset+7] - data[i+4 .. i+7] (calculate difference)
        EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[next 4], q4 = q4 * w (multiply difference by twiddle)
        EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[next 4], q4 = q4 * w (multiply difference by next twiddle)
        ee.vmin.s16.ST.INCP q4, a15,  q6, q3, q6  // data[next 4 offset] = q4 (store difference), update minimum in q6 using q3 values
        ee.vmax.s16  q7, q3, q7  // Update maximum in q7 using q3 values
        ee.vmin.s16  q6, q4, q6  // Update minimum in q6 using q4 values
        ee.vmax.s16  q7, q4, q7  // Update maximum in q7 using q4 values
   .loop1_end_stageN:
    add.n a10,a10,a12   // Move to the next group of butterflies (a10 += offset - 16)
    bne a6, a11, .loop0_start_stageN  // Continue outer loop if we haven't processed all groups

    // Find global minimum and maximum across all vector elements in q6 and q7
    ee.vunzip.16    q6, q4    // Unzip q6 and q4 to prepare for reduction
    ee.vmin.s16     q6, q4, q6 // minimum x4 (reduce within pairs)
    ee.vunzip.16    q6, q4
    ee.vmin.s16     q6, q4, q6 // minimum x2 (reduce further)
    ee.vunzip.16    q6, q4
    ee.vmin.s16     q6, q4, q6 // minimum x1 (reduce to single minimum in q6)
    ee.vunzip.16    q7, q4
    ee.vmax.s16     q7, q4, q7 // maximum x4 (reduce within pairs)
    ee.vunzip.16    q7, q4
    ee.vmax.s16     q7, q4, q7 // maximum x2 (reduce further)
    ee.vunzip.16    q7, q4
    ee.vmax.s16     q7, q4, q7 // maximum x1 (reduce to single maximum in q7)
    ee.movi.32.a	q6, a5, 0 // Move the minimum value from q6 (low 32 bits) to register a5
    slli 	        a5, a5, 16  // Shift minimum left by 16 bits (likely to get the high 16 bits of the int16_t value)
    abs 	        a5, a5      // Take absolute value of the minimum
    ee.movi.32.a	q7, a14, 0 // Move the maximum value from q7 (low 32 bits) to register a14
    slli 	        a14, a14, 16  // Shift maximum left by 16 bits
    max             a5, a5, a14 // a5 = max(abs(minimum), maximum) (find the largest absolute value)

    beqz    a5, .end  // If the largest absolute value is zero, no scaling needed, jump to end
    l32r	a15, .LC1_MX_32  // Load the maximum 32-bit signed integer value (0x7FFFFFFF)
    quos    a5, a15, a5  // a5 = 0x7FFFFFFF / (max_abs_value << 16) (calculate scaling factor)
    bnei    a5,1,.normalize_stageN_0  // If scaling factor is not 1, jump to normalization
    // Scaling factor is 1, use 0.5 for normalization (no actual scaling, just normalization by 2)
    l32r a14,.Q15_05
    ee.movi.32.q	q5, a14, 0  // Load 0.5 into q5
    ee.movi.32.q	q5, a14, 1 
    ee.movi.32.q	q5, a14, 2 
    ee.movi.32.q	q5, a14, 3
    j .normalize_stageN  // Jump to normalization loop
.normalize_stageN_0: 
    blti a5,4,.end_normalize_stageN  // If scaling factor is less than 4, skip normalization (values are small enough)
    slli            a15, a5, 16    // a15 = scaling_factor << 16
    or              a15, a5, a15   // a15 = scaling_factor | (scaling_factor << 16) (create 16-bit repeating pattern)
    ee.movi.32.q	q5, a15, 0  // Load the scaling factor into all slots of q5
    ee.movi.32.q	q5, a15, 1 
    ee.movi.32.q	q5, a15, 2 
    ee.movi.32.q	q5, a15, 3
    ssai            1         // Set SAR register to 1 (for the normalization multiplication shift)
.normalize_stageN:
    mov.n a10,a2    // a10 = &data[0] (load start of data array)
    mov.n a11,a2    // a11 = &data[0] (load start of data array for output)
    ee.vld.128.ip   q1, a10, 16  // Load first 128 bits of data into q1, increment a10
    loopnez a7, .stageN_end_normalize  // Loop a7 times (size/4 iterations)
        ee.vmul.s16.ld.incp     q1, a10, q2, q1, q5  // Multiply q1 by scaling factor q5, store result in q2, load next data into q1, increment a10
        ee.vst.128.ip   q2, a11, 16  // Store the scaled result from q2, increment a11
    .stageN_end_normalize:
    ssai            15         // Set SAR register back to 15
.end_normalize_stageN:

    slli a8,a8,1     // s *= 2 (double the stage size for the next iteration)
    mov.n a3,a9      // Update twiddle factor pointer (a9 was advanced during the loop)
    bnei a4,32,.stageN  // Continue stages N if offset is not 32 (i.e., not the 3rd stage from the end)

.stage3:             // Stage where butterfly size is 4 (special case)
    ee.zero.q q6 // Initialize minimum value accumulator q6 to zero
    ee.zero.q q7 // Initialize maximum value accumulator q7 to zero

    mov.n a10,a2        // a10 = &data[0] (start of data for this stage)
    addi.n a11,a10,16   // a11 = &data[offset] (pointer to second operand, offset=4*4=16 bytes for 4 complex numbers)
    mov.n a14,a10       // a14 = &data[i] (save pointer for writing result back)
    mov.n a15,a11       // a15 = &data[i+offset] (save pointer for writing result back)

    EE.VLD.128.IP q1, a10, 32                   // q1 = data[i .. i+3] (load first 4 complex numbers), increment a10 by 32 (skip next 4 for now)
    EE.VLD.128.IP q2, a11, 32                   // q2 = data[i+offset .. i+offset+3] (load second 4 complex numbers), increment a11 by 32
    EE.VADDS.S16.LD.INCP q0, a9, q3, q1, q2     // q0 = w[twiddle_index .. twiddle_index+3], q3 = data[i .. i+3] + data[i+offset .. i+offset+3]
    EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1    // data[i .. i+3] = q3 (store sum), q4 = data[i+offset .. i+offset+3] - data[i .. i+3] (calculate difference)
    EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[i+8 .. i+11] (reload from new address), q4 = q4 * w (multiply difference by twiddle)
    EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[i+offset+8 .. i+offset+11] (reload from new address), q4 = q4 * w (multiply difference by next twiddle)
    ee.vmin.s16.ST.INCP q4, a15,  q6, q3, q6  // data[next 4 offset] = q4 (store difference), update minimum in q6 using q3 values
    ee.vmax.s16.LD.INCP q5, a10,    q7, q3, q7  // Update maximum in q7 using q3 values, load next data into q5
    ee.vmin.s16.LD.INCP q5, a11,    q6, q4, q6  // Update minimum in q6 using q4 values, load next data into q5
    ee.vmax.s16.LD.INCP q5, a14,    q7, q4, q7  // Update maximum in q7 using q4 values, load next data into q5
 
    addi.n a5,a8,-1     // Calculate number of groups for this stage (a8 butterflies per group, so a8-1 iterations)
    loopnez a5, .loop1_end_stage3  // Inner loop: process remaining groups for this stage
        EE.VADDS.S16.LD.INCP q5, a15, q3, q1, q2    // a14+=16, q3 = data[next 4] + data[next 4] (from updated a14)
        EE.VSUBS.S16.ST.INCP q3, a14, q4, q2, q1    // data[next 4] = q3 (store sum), q4 = data[next 4 offset] - data[next 4] (calculate difference)
        EE.CMUL.S16.LD.INCP q1, a10, q4, q4, q0, 0  // q1 = data[next 4], q4 = q4 * w (multiply difference by twiddle)
        EE.CMUL.S16.LD.INCP q2, a11, q4, q4, q0, 1  // q2 = data[next 4], q4 = q4 * w (multiply difference by next twiddle)
        
        ee.vmin.s16.ST.INCP q4, a15,  q6, q3, q6  // data[next 4 offset] = q4 (store difference), update minimum in q6 using q3 values
        ee.vmax.s16.LD.INCP q5, a10,    q7, q3, q7  // Update maximum in q7 using q3 values, load next data into q5
        ee.vmin.s16.LD.INCP q5, a11,    q6, q4, q6  // Update minimum in q6 using q4 values, load next data into q5
        ee.vmax.s16.LD.INCP q5, a14,    q7, q4, q7  // Update maximum in q7 using q4 values, load next data into q5
    .loop1_end_stage3:

    // Find global minimum and maximum across all vector elements in q6 and q7 (similar to stageN)
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

    beqz    a5, .end  // If the largest absolute value is zero, no scaling needed, jump to end
    l32r	a15, .LC1_MX_32  // Load the maximum 32-bit signed integer value (0x7FFFFFFF)
    quos    a5, a15, a5  // a5 = 0x7FFFFFFF / (max_abs_value << 16) (calculate scaling factor)
    bnei    a5,1,.normalize_stage2_0  // If scaling factor is not 1, jump to normalization
    // Scaling factor is 1, use 0.5 for normalization (no actual scaling, just normalization by 2)
    l32r a14,.Q15_05
    ee.movi.32.q	q5, a14, 0  // Load 0.5 into q5
    ee.movi.32.q	q5, a14, 1 
    ee.movi.32.q	q5, a14, 2 
    ee.movi.32.q	q5, a14, 3
    j .normalize_stage2  // Jump to normalization loop
.normalize_stage2_0: 
    blti a5,4,.stage2  // If scaling factor is less than 4, skip normalization for this stage and go to stage2
    slli            a15, a5, 16    // a15 = scaling_factor << 16
    or              a15, a5, a15   // a15 = scaling_factor | (scaling_factor << 16) (create 16-bit repeating pattern)
    ee.movi.32.q	q5, a15, 0  // Load the scaling factor into all slots of q5
    ee.movi.32.q	q5, a15, 1 
    ee.movi.32.q	q5, a15, 2 
    ee.movi.32.q	q5, a15, 3
    ssai            1         // Set SAR register to 1 (for the normalization multiplication shift)
.normalize_stage2:
    mov.n a10,a2    // a10 = &data[0] (load start of data array)
    mov.n a11,a2    // a11 = &data[0] (load start of data array for output)
    ee.vld.128.ip   q1, a10, 16  // Load first 128 bits of data into q1, increment a10
    loopnez a7, .stage2_end_normalize  // Loop a7 times (size/4 iterations)
        ee.vmul.s16.ld.incp     q1, a10, q2, q1, q5  // Multiply q1 by scaling factor q5, store result in q2, load next data into q1, increment a10
        ee.vst.128.ip   q2, a11, 16  // Store the scaled result from q2, increment a11
    .stage2_end_normalize:
    ssai            15         // Set SAR register back to 15 
////////////////////////////////////////////////////////////////////////// 
 .stage2:            // Stage where butterfly size is 2 (special case)
    ee.zero.q q6 // Initialize minimum value accumulator q6 to zero
    ee.zero.q q7 // Initialize maximum value accumulator q7 to zero

    mov.n a10,a2    // a10 = &data[0] (start of data for this stage)
    mov.n a11,a2    // a11 = &data[0] (start of data for output in this stage)

    EE.VLD.H.64.IP q0, a9, 0                    // q0hi = w[0 1] (load first two twiddle factors into high part of q0)
    EE.VLD.128.IP q1, a10, 16                   // q1 = data[i i+1 i+offset i+offset+1] (load 2 complex pairs)
    EE.VLD.128.IP q2, a10, 16     //  q2 = data[i+2 i+3 i+offset+2 i+offset+3] (load next 2 complex pairs)
    loopnez a8, .loop1_end_stage2  // Loop a8 times (number of groups in this stage)
        EE.FFT.R2BF.S16 q3, q4, q1, q2, 0  // Perform radix-2 butterfly operation: q3 = out1, q4 = out2
        EE.CMUL.S16.LD.INCP q1, a10, q3, q3, q0, 1  // q1 = data[next pair], q3hi = q3hi * w[0 1] (multiply high parts of result by twiddle factors)
        EE.CMUL.S16.LD.INCP q2, a10, q4, q4, q0, 1 // q2 = data[next pair], q4hi = q4hi * w[0 1] (multiply high parts of result by twiddle factors)
        ee.vmin.s16.ST.INCP q3, a11, q6, q3, q6  // data[next pair] = q3 (store result), update minimum in q6 using q3 values
        ee.vmax.s16.ST.INCP q4, a11, q7, q3, q7  // data[next pair] = q4 (store result), update maximum in q7 using q3 values
        ee.vmin.s16  q6, q4, q6  // Update minimum in q6 using q4 values
        ee.vmax.s16  q7, q4, q7  // Update maximum in q7 using q4 values
    .loop1_end_stage2:
////////////////////////////////////
    // Find global minimum and maximum across all vector elements in q6 and q7 (similar to stageN/stage3)
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
    //S32I.N a5, a2, 0  // (Commented out: Store max absolute value to the beginning of data array)
    //retw.n             // (Commented out: Return immediately)


    beqz    a5, .end  // If the largest absolute value is zero, no scaling needed, jump to end
    l32r	a15, .LC1_MX_32  // Load the maximum 32-bit signed integer value (0x7FFFFFFF)
    quos    a5, a15, a5  // a5 = 0x7FFFFFFF / (max_abs_value << 16) (calculate scaling factor)
    bnei    a5,1,.normalize_stage1_0  // If scaling factor is not 1, jump to normalization
    // Scaling factor is 1, use 0.5 for normalization (no actual scaling, just normalization by 2)
    l32r a14,.Q15_05
    ee.movi.32.q	q5, a14, 0  // Load 0.5 into q5
    ee.movi.32.q	q5, a14, 1 
    ee.movi.32.q	q5, a14, 2 
    ee.movi.32.q	q5, a14, 3
    j .normalize_stage1  // Jump to normalization loop
.normalize_stage1_0: 
    blti a5,4,.stage1  // If scaling factor is less than 4, skip normalization for this stage and go to stage1
    slli            a15, a5, 16    // a15 = scaling_factor << 16
    or              a15, a5, a15   // a15 = scaling_factor | (scaling_factor << 16) (create 16-bit repeating pattern)
    ee.movi.32.q	q5, a15, 0  // Load the scaling factor into all slots of q5
    ee.movi.32.q	q5, a15, 1 
    ee.movi.32.q	q5, a15, 2 
    ee.movi.32.q	q5, a15, 3
    ssai            1         // Set SAR register to 1 (for the normalization multiplication shift)
.normalize_stage1:
    mov.n a10,a2    // a10 = &data[0] (load start of data array)
    mov.n a11,a2    // a11 = &data[0] (load start of data array for output)
    ee.vld.128.ip   q1, a10, 16  // Load first 128 bits of data into q1, increment a10
    loopnez a7, .stage1_end_normalize  // Loop a7 times (size/4 iterations)
        ee.vmul.s16.ld.incp     q1, a10, q2, q1, q5  // Multiply q1 by scaling factor q5, store result in q2, load next data into q1, increment a10
        ee.vst.128.ip   q2, a11, 16  // Store the scaled result from q2, increment a11
    .stage1_end_normalize:
    ssai            15         // Set SAR register back to 15 
////////////////////////////////////
.stage1:            // Stage where butterfly size is 1 (final stage)
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

.end:
    retw.n  // Return from function (windowed return)
