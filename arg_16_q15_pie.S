    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access, often required for time-critical code like ISRs or frequently called functions)
    .type   arg_16_q15_pie,@function
    .align  4
    .global arg_16_q15_pie
arg_16_q15_pie: 
// in   - a2  (pointer to input complex vector)
// out  - a3  (pointer to output argument vector)
// size - a4  (number of complex numbers to process)
// tan  - a5  (pointer to lookup table of tangent values)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI
    
    MOVI.N a9, 13     // Load loop counter for CORDIC iterations (a9 = 13)
    MOVI.N a7, 0      // Load constant 0 (a7 = 0, often used as zero vector or immediate)
    MOVI.N a6, -1     // Load constant -1 (a6 = -1, often used as mask)
    EE.ZERO.Q q7      // Clear vector register q7 (used as zero vector)
    MOV.N a10, a5     // Save the original tan pointer (a5) to a10

    slli 	a4, a4, 1           // size * 2 (calculate byte offset for end of output array)
    add.n 	a4, a4, a3          // a4 = out + size*2 (calculate end address of output array)

    EE.VLD.128.IP q0, a2, 16  // Load 128 bits (4 complex_q15 = 8 int16_t) from input into q0, increment a2 by 16
    EE.VLD.128.IP q1, a2, 16  // Load next 128 bits (4 complex_q15 = 8 int16_t) from input into q1, increment a2 by 16

    .big_loop_start:          // Outer loop: process 8 complex numbers per iteration
        EE.VUNZIP.16 q0, q1         // Interleave q0 and q1 to separate real (x) into q0, imaginary (y) into q1
        EE.ZERO.Q q6            // Initialize sum accumulator q6 to zero (for CORDIC angle calculation)
        EE.MOV.S16.QACC q0      // Move q0 (real parts) to Q accumulator (for later sign extraction)
        EE.VRELU.S16 q0, a6, a7     // Take absolute value of q0 (real parts), using mask a6 (-1) and immediate a7 (0) -> q0 now holds |x|
        MOVI.N a8, 0            // Initialize inner loop counter a8 to 0
        loopnez a9, .loop_cordic_end  // Inner loop: perform 13 CORDIC iterations
            // step 0
            EE.VCMP.LT.S16 q5, q1, q7  // Compare q1 (y) < 0 (q7=0), result mask in q5 (for negative y)
            EE.VCMP.GT.S16 q4, q1, q7  // Compare q1 (y) > 0 (q7=0), result mask in q4 (for positive y)
            // xn = x1 + (y >> a8) or xn = x1 - (y >> a8)  (depending on sign of y)
            WSR.SAR a8              // Set shift amount register for vector shifts
            EE.VMUL.S16 q3, q1, q4 // q3 = y * (mask for positive y)
            EE.VSUBS.S16 q2, q0, q3 // q2 = |x| - (y if y>0 else 0)  (subtract if y positive)
            EE.VMUL.S16 q3, q1, q5 // q3 = y * (mask for negative y)
            EE.VADDS.S16 q2, q2, q3 // q2 = q2 + (y if y<0 else 0)  (add if y negative) -> xn
            // y = y - (x1 >> a8) or y = y + (x1 >> a8) (depending on sign of x)
            EE.VMUL.S16 q3, q0, q4 // q3 = |x| * (mask for positive y)  (use sign of original y)
            EE.VADDS.S16 q1, q1, q3 // q1 = y + (|x| if y>0 else 0)  (add if original y was positive)
            EE.VMUL.S16 q3, q0, q5 // q3 = |x| * (mask for negative y)
            EE.VSUBS.S16 q1, q1, q3 // q1 = q1 - (|x| if y<0 else 0)  (subtract if original y was negative)
            // x1 = xn (update x register)
            EE.VADDS.S16 q0, q2, q7     // q0 = xn (new x value)
            // sum += tan[a8] or sum -= tan[a8] (update angle accumulator based on rotation direction)
            WSR.SAR a7              // Set shift amount to 0 (or potentially reset SAR register)
            EE.VLDBC.16.IP q2, a5, 2    // Load 8 copies of tan[a8] into q2, increment a5 by 2 (points to next tan value)
            EE.VMUL.S16 q4, q2, q4      // q4 = tan[a8] * (mask for positive y)
            EE.VSUBS.S16 q6, q6, q4     // sum = sum - (tan[a8] if y>0 else 0)  (subtract if y positive)
            EE.VMUL.S16 q5, q2, q5      // q5 = tan[a8] * (mask for negative y)
            EE.VADDS.S16 q6, q6, q5     // sum = sum + (tan[a8] if y<0 else 0)  (add if y negative)

            ADDI.N a8,a8,1          // Increment inner loop counter a8
        .loop_cordic_end:           // End of inner CORDIC loop
        
        EE.SRCMB.S16.QACC q3, a7, 0 // Extract sign of original x from Q accumulator -> q3 (0 or -1)
        EE.VLDBC.16.IP q2, a5, 2    // Load 8 copies of -16383 (PI) into q2, increment a5
        EE.VLDBC.16.IP q1, a5, 2    // Load 8 copies of 1 into q1, increment a5
        EE.VCMP.LT.S16 q5, q3, q7   // Check if original x was negative (q3 < 0) -> q5 (mask for neg x)
        EE.VCMP.LT.S16 q0, q6, q7   // Check if CORDIC sum was negative (q6 < 0) -> q0 (mask for neg sum)
        EE.VMUL.S16 q4, q5, q0      // Combine masks: true if both x and sum were negative -> q4
        EE.VMUL.S16 q4, q2, q4      // If both neg, load -16383 into q4 (correction for quadrant)

        EE.VADDS.S16 q0, q0, q1     // Invert mask q0 (not neg sum) -> q0
        EE.VMUL.S16 q0, q5, q0      // Combine: true if x neg and sum not neg -> q0
        EE.VMUL.S16 q0, q2, q0      // If x neg and sum pos, load 16383 into q0 (correction for quadrant)
        EE.VADDS.S16 q3, q0, q4     // Combine the two possible corrections (-16383 or +16383) -> q3

        EE.VCMP.EQ.S16 q4, q3, q7   // Check if q3 (correction) is zero (i.e., original x was non-negative) -> q4 (mask for not neg x)
        EE.VLDBC.16 q5, a5          // Load 8 copies of 2 into q5 (used for sign adjustment)
        EE.VMUL.S16.LD.INCP q0, a2, q4, q4, q5  // Load next 8 real parts from input, multiply by q4 (0 if x>=0, else 2), add to q4, increment a2 -> q0 (this seems complex, likely part of sign handling or input fetch)
        EE.VADDS.S16 q4, q4, q1     // Add 1 to q4 (likely creating a sign mask based on x's sign)
        EE.VMUL.S16.LD.INCP q1, a2, q6, q6, q4  // Load next 8 imag parts from input, multiply by q6 (sum), add to q4 (sign), increment a2 -> q1 (this also seems complex, likely part of input fetch or intermediate calc)

        EE.VSUBS.S16 q6, q3, q6     // Final angle: correction (q3) - CORDIC sum (q6) -> q6

        EE.VST.128.IP q6, a3, 16 // Save calculated argument (q6) to output array, increment a3 by 16
        MOV.N a5, a10             // Restore the original tan pointer for the next outer loop iteration

        BLT a3, a4, .big_loop_start  // Branch back to .big_loop_start if a3 (output pointer) is less than a4 (end address)

    retw.n  // Return from function (windowed return)
