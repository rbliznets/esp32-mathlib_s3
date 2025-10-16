    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access)
    .text
    .align  4
	.literal_position
	.literal	.PLUS_PI, 16383      // Literal for +Pi in Q15 format
	.literal	.MINUS_PI, -16383    // Literal for -Pi in Q15 format
	.literal	.PLUS_PI2, 8191      // Literal for +Pi/2 in Q15 format
	.literal	.MINUS_PI2, -8191    // Literal for -Pi/2 in Q15 format
    .literal	.COS_FACTOR, 19894   // Initial gain factor for cosine calculation in CORDIC
	.literal	.SAT_MAX_16, 32767   // Literal for maximum value of int16_t
	.literal	.SAT_MIN_16, -32768  // Literal for minimum value of int16_t
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   atan2_q15_s3,@function
    .align  4
    .global atan2_q15_s3
atan2_q15_s3: 
// y   - a2  (input y coordinate)
// x   - a3  (input x coordinate)
// &tan - a4 (pointer to lookup table of tangent values for CORDIC)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    abs 	a5, a3       // x1 = absolute value of x (a5 = abs(a3))
    movi.n  a6, 0        // sum = 0 (initialize angle accumulator a6)

    // 0th CORDIC iteration
    l16si   a8, a4, 0   // Load tan[0] from lookup table into a8
    bgez a2, .atan2_step0  // Branch if y >= 0
        // y < 0 (adjust rotation direction)
        sub a7,a5,a2    // xn = x1 - y (CORDIC x update)
        add.n a2,a2,a5  // y = y + x1 (CORDIC y update)
        sub a6,a6,a8  // sum -= tan[0] (adjust angle accumulator)
        j .atan2_step0_end  // Jump to end of step 0
    .atan2_step0:
        // y >= 0 (adjust rotation direction)
        add.n a7,a5,a2  // xn = x1 + y (CORDIC x update)
        sub a2,a2,a5    // y = y - x1 (CORDIC y update)
        add.n a6,a6,a8  // sum += tan[0] (adjust angle accumulator)
    .atan2_step0_end:
        beqz a2,.atan2_step12_end  // Exit early if y is now zero
        mov.n a5,a7     // x1 = xn (update x for next iteration)

    // 1st CORDIC iteration
    l16si a8, a4, 2     // Load tan[1] from lookup table into a8
    srai a9, a5, 1      // a9 = x1 >> 1 (pre-shifted x for CORDIC)
    srai a10, a2, 1     // a10 = y >> 1 (pre-shifted y for CORDIC)
    bgez a2, .atan2_step1  // Branch if y >= 0
        // y < 0 (adjust rotation direction)
        sub a7,a5,a10    // xn = x1 - (y >> 1) (CORDIC x update)
        add.n a2,a2,a9  // y = y + (x1 >> 1) (CORDIC y update)
        sub a6,a6,a8  // sum -= tan[1] (adjust angle accumulator)
        j .atan2_step1_end  // Jump to end of step 1
    .atan2_step1:
        // y >= 0 (adjust rotation direction)
        add.n a7,a5,a10  // xn = x1 + (y >> 1) (CORDIC x update)
        sub a2,a2,a9    // y = y - (x1 >> 1) (CORDIC y update)
        add.n a6,a6,a8  // sum += tan[1] (adjust angle accumulator)
    .atan2_step1_end:
        beqz a2,.atan2_step12_end  // Exit early if y is now zero
        mov.n a5,a7    // x1 = xn (update x for next iteration)

    // 2nd CORDIC iteration
    l16si a8, a4, 4     // Load tan[2] from lookup table into a8
    srai a9, a5, 2      // a9 = x1 >> 2 (pre-shifted x for CORDIC)
    srai a10, a2, 2     // a10 = y >> 2 (pre-shifted y for CORDIC)
    bgez a2, .atan2_step2  // Branch if y >= 0
        // y < 0 (adjust rotation direction)
        sub a7,a5,a10    // xn = x1 - (y >> 2) (CORDIC x update)
        add.n a2,a2,a9  // y = y + (x1 >> 2) (CORDIC y update)
        sub a6,a6,a8  // sum -= tan[2] (adjust angle accumulator)
        j .atan2_step2_end  // Jump to end of step 2
    .atan2_step2:
        // y >= 0 (adjust rotation direction)
        add.n a7,a5,a10  // xn = x1 + (y >> 2) (CORDIC x update)
        sub a2,a2,a9    // y = y - (x1 >> 2) (CORDIC y update)
        add.n a6,a6,a8  // sum += tan[2] (adjust angle accumulator)
    .atan2_step2_end:
        beqz a2,.atan2_step12_end  // Exit early if y is now zero
        mov.n a5,a7    // x1 = xn (update x for next iteration)

    // 3rd CORDIC iteration
    l16si a8, a4, 6     // Load tan[3] from lookup table into a8
    srai a9, a5, 3      // a9 = x1 >> 3 (pre-shifted x for CORDIC)
    srai a10, a2, 3     // a10 = y >> 3 (pre-shifted y for CORDIC)
    bgez a2, .atan2_step3  // Branch if y >= 0
        // y < 0 (adjust rotation direction)
        sub a7,a5,a10    // xn = x1 - (y >> 3) (CORDIC x update)
        add.n a2,a2,a9  // y = y + (x1 >> 3) (CORDIC y update)
        sub a6,a6,a8  // sum -= tan[3] (adjust angle accumulator)
        j .atan2_step3_end  // Jump to end of step 3
    .atan2_step3:
        // y >= 0 (adjust rotation direction)
        add.n a7,a5,a10  // xn = x1 + (y >> 3) (CORDIC x update)
        sub a2,a2,a9    // y = y - (x1 >> 3) (CORDIC y update)
        add.n a6,a6,a8  // sum += tan[3] (adjust angle accumulator)
    .atan2_step3_end:
        beqz a2,.atan2_step12_end  // Exit early if y is now zero
        mov.n a5,a7    // x1 = xn (update x for next iteration)

    // 4th CORDIC iteration
    l16si a8, a4, 8     // Load tan[4] from lookup table into a8
    srai a9, a5, 4      // a9 = x1 >> 4 (pre-shifted x for CORDIC)
    srai a10, a2, 4     // a10 = y >> 4 (pre-shifted y for CORDIC)
    bgez a2, .atan2_step4  // Branch if y >= 0
        // y < 0 (adjust rotation direction)
        sub a7,a5,a10    // xn = x1 - (y >> 4) (CORDIC x update)
        add.n a2,a2,a9  // y = y + (x1 >> 4) (CORDIC y update)
        sub a6,a6,a8  // sum -= tan[4] (adjust angle accumulator)
        j .atan2_step4_end  // Jump to end of step 4
    .atan2_step4:
        // y >= 0 (adjust rotation direction)
        add.n a7,a5,a10  // xn = x1 + (y >> 4) (CORDIC x update)
        sub a2,a2,a9    // y = y - (x1 >> 4) (CORDIC y update)
        add.n a6,a6,a8  // sum += tan[4] (adjust angle accumulator)
    .atan2_step4_end:
        beqz a2,.atan2_step12_end  // Exit early if y is now zero
        mov.n a5,a7    // x1 = xn (update x for next iteration)

    // 5th CORDIC iteration
    l16si a8, a4, 10     // Load tan[5] from lookup table into a8
    srai a9, a5, 5      // a9 = x1 >> 5 (pre-shifted x for CORDIC)
    srai a10, a2, 5     // a10 = y >> 5 (pre-shifted y for CORDIC)
    bgez a2, .atan2_step5  // Branch if y >= 0
        // y < 0 (adjust rotation direction)
        sub a7,a5,a10    // xn = x1 - (y >> 5) (CORDIC x update)
        add.n a2,a2,a9  // y = y + (x1 >> 5) (CORDIC y update)
        sub a6,a6,a8  // sum -= tan[5] (adjust angle accumulator)
        j .atan2_step5_end  // Jump to end of step 5
    .atan2_step5:
        // y >= 0 (adjust rotation direction)
        add.n a7,a5,a10  // xn = x1 + (y >> 5) (CORDIC x update)
        sub a2,a2,a9    // y = y - (x1 >> 5) (CORDIC y update)
        add.n a6,a6,a8  // sum += tan[5] (adjust angle accumulator)
    .atan2_step5_end:
        beqz a2,.atan2_step12_end  // Exit early if y is now zero
        mov.n a5,a7    // x1 = xn (update x for next iteration)

    // 6th CORDIC iteration
    l16si a8, a4, 12     // Load tan[6] from lookup table into a8
    srai a9, a5, 6      // a9 = x1 >> 6 (pre-shifted x for CORDIC)
    srai a10, a2, 6     // a10 = y >> 6 (pre-shifted y for CORDIC)
    bgez a2, .atan2_step6  // Branch if y >= 0
        // y < 0 (adjust rotation direction)
        sub a7,a5,a10    // xn = x1 - (y >> 6) (CORDIC x update)
        add.n a2,a2,a9  // y = y + (x1 >> 6) (CORDIC y update)
        sub a6,a6,a8  // sum -= tan[6] (adjust angle accumulator)
        j .atan2_step6_end  // Jump to end of step 6
    .atan2_step6:
        // y >= 0 (adjust rotation direction)
        add.n a7,a5,a10  // xn = x1 + (y >> 6) (CORDIC x update)
        sub a2,a2,a9    // y = y - (x1 >> 6) (CORDIC y update)
        add.n a6,a6,a8  // sum += tan[6] (adjust angle accumulator)
    .atan2_step6_end:
        beqz a2,.atan2_step12_end  // Exit early if y is now zero
        mov.n a5,a7    // x1 = xn (update x for next iteration)

    // 7th CORDIC iteration
    l16si a8, a4, 14     // Load tan[7] from lookup table into a8
    srai a9, a5, 7      // a9 = x1 >> 7 (pre-shifted x for CORDIC)
    srai a10, a2, 7     // a10 = y >> 7 (pre-shifted y for CORDIC)
    bgez a2, .atan2_step7  // Branch if y >= 0
        // y < 0 (adjust rotation direction)
        sub a7,a5,a10    // xn = x1 - (y >> 7) (CORDIC x update)
        add.n a2,a2,a9  // y = y + (x1 >> 7) (CORDIC y update)
        sub a6,a6,a8  // sum -= tan[7] (adjust angle accumulator)
        j .atan2_step7_end  // Jump to end of step 7
    .atan2_step7:
        // y >= 0 (adjust rotation direction)
        add.n a7,a5,a10  // xn = x1 + (y >> 7) (CORDIC x update)
        sub a2,a2,a9    // y = y - (x1 >> 7) (CORDIC y update)
        add.n a6,a6,a8  // sum += tan[7] (adjust angle accumulator)
    .atan2_step7_end:
        beqz a2,.atan2_step12_end  // Exit early if y is now zero
        mov.n a5,a7    // x1 = xn (update x for next iteration)

    // 8th CORDIC iteration
    l16si a8, a4, 16     // Load tan[8] from lookup table into a8
    srai a9, a5, 8      // a9 = x1 >> 8 (pre-shifted x for CORDIC)
    srai a10, a2, 8     // a10 = y >> 8 (pre-shifted y for CORDIC)
    bgez a2, .atan2_step8  // Branch if y >= 0
        // y < 0 (adjust rotation direction)
        sub a7,a5,a10    // xn = x1 - (y >> 8) (CORDIC x update)
        add.n a2,a2,a9  // y = y + (x1 >> 8) (CORDIC y update)
        sub a6,a6,a8  // sum -= tan[8] (adjust angle accumulator)
        j .atan2_step8_end  // Jump to end of step 8
    .atan2_step8:
        // y >= 0 (adjust rotation direction)
        add.n a7,a5,a10  // xn = x1 + (y >> 8) (CORDIC x update)
        sub a2,a2,a9    // y = y - (x1 >> 8) (CORDIC y update)
        add.n a6,a6,a8  // sum += tan[8] (adjust angle accumulator)
    .atan2_step8_end:
        beqz a2,.atan2_step12_end  // Exit early if y is now zero
        mov.n a5,a7    // x1 = xn (update x for next iteration)

    // 9th CORDIC iteration
    l16si a8, a4, 18     // Load tan[9] from lookup table into a8
    srai a9, a5, 9      // a9 = x1 >> 9 (pre-shifted x for CORDIC)
    srai a10, a2, 9     // a10 = y >> 9 (pre-shifted y for CORDIC)
    bgez a2, .atan2_step9  // Branch if y >= 0
        // y < 0 (adjust rotation direction)
        sub a7,a5,a10    // xn = x1 - (y >> 9) (CORDIC x update)
        add.n a2,a2,a9  // y = y + (x1 >> 9) (CORDIC y update)
        sub a6,a6,a8  // sum -= tan[9] (adjust angle accumulator)
        j .atan2_step9_end  // Jump to end of step 9
    .atan2_step9:
        // y >= 0 (adjust rotation direction)
        add.n a7,a5,a10  // xn = x1 + (y >> 9) (CORDIC x update)
        sub a2,a2,a9    // y = y - (x1 >> 9) (CORDIC y update)
        add.n a6,a6,a8  // sum += tan[9] (adjust angle accumulator)
    .atan2_step9_end:
        beqz a2,.atan2_step12_end  // Exit early if y is now zero
        mov.n a5,a7    // x1 = xn (update x for next iteration)

    // 10th CORDIC iteration
    l16si a8, a4, 20     // Load tan[10] from lookup table into a8
    srai a9, a5, 10      // a9 = x1 >> 10 (pre-shifted x for CORDIC)
    srai a10, a2, 10     // a10 = y >> 10 (pre-shifted y for CORDIC)
    bgez a2, .atan2_step10  // Branch if y >= 0
        // y < 0 (adjust rotation direction)
        sub a7,a5,a10    // xn = x1 - (y >> 10) (CORDIC x update)
        add.n a2,a2,a9  // y = y + (x1 >> 10) (CORDIC y update)
        sub a6,a6,a8  // sum -= tan[10] (adjust angle accumulator)
        j .atan2_step10_end  // Jump to end of step 10
    .atan2_step10:
        // y >= 0 (adjust rotation direction)
        add.n a7,a5,a10  // xn = x1 + (y >> 10) (CORDIC x update)
        sub a2,a2,a9    // y = y - (x1 >> 10) (CORDIC y update)
        add.n a6,a6,a8  // sum += tan[10] (adjust angle accumulator)
    .atan2_step10_end:
        beqz a2,.atan2_step12_end  // Exit early if y is now zero
        mov.n a5,a7    // x1 = xn (update x for next iteration)

    // 11th CORDIC iteration
    l16si a8, a4, 22     // Load tan[11] from lookup table into a8
    srai a9, a5, 11      // a9 = x1 >> 11 (pre-shifted x for CORDIC)
    srai a10, a2, 11     // a10 = y >> 11 (pre-shifted y for CORDIC)
    bgez a2, .atan2_step11  // Branch if y >= 0
        // y < 0 (adjust rotation direction)
        sub a7,a5,a10    // xn = x1 - (y >> 11) (CORDIC x update)
        add.n a2,a2,a9  // y = y + (x1 >> 11) (CORDIC y update)
        sub a6,a6,a8  // sum -= tan[11] (adjust angle accumulator)
        j .atan2_step11_end  // Jump to end of step 11
    .atan2_step11:
        // y >= 0 (adjust rotation direction)
        add.n a7,a5,a10  // xn = x1 + (y >> 11) (CORDIC x update)
        sub a2,a2,a9    // y = y - (x1 >> 11) (CORDIC y update)
        add.n a6,a6,a8  // sum += tan[11] (adjust angle accumulator)
    .atan2_step11_end:
        beqz a2,.atan2_step12_end  // Exit early if y is now zero

    // 12th CORDIC iteration (last one, only updates angle accumulator)
    l16si a8, a4, 24     // Load tan[12] from lookup table into a8
    bgez a2, .atan2_step12  // Branch if y >= 0
        // y < 0 (adjust rotation direction)
        sub a6,a6,a8  // sum -= tan[12] (adjust angle accumulator)
        j .atan2_step12_end  // Jump to end of step 12
    .atan2_step12:
        // y >= 0 (adjust rotation direction)
        add.n a6,a6,a8  // sum += tan[12] (adjust angle accumulator)
    .atan2_step12_end:
   
    bgez a3,.end_atan2_q15_s3  // Branch if x >= 0 (no quadrant correction needed)
        // x < 0 (apply quadrant correction)
        bgez a6,.atan2_last_check  // Branch if sum >= 0
            // sum < 0 (result is in the 3rd quadrant, angle = PI - calculated_angle)
            l32r a7,.MINUS_PI  // Load -PI constant
            j .end1_atan2_q15_s3  // Jump to final calculation
        .atan2_last_check:
            // sum >= 0 (result is in the 2nd quadrant, angle = PI - calculated_angle)
            l32r a7,.PLUS_PI  // Load +PI constant
        .end1_atan2_q15_s3:
        sub a2,a7,a6  // Calculate final result: PI - sum or -PI - sum
        retw.n        // Return from function (windowed return)
    .end_atan2_q15_s3:  // x >= 0 case
    mov.n a2,a6  // Result is the calculated sum
    retw.n       // Return from function (windowed return)
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   sincos_q15_s3,@function
    .align  4
    .global sincos_q15_s3
sincos_q15_s3: 
// angle   - a2  (input angle in Q15 format)
// &cos   - a3   (pointer to output cosine value)
// &sin   - a4   (pointer to output sine value)
// &tan   - a5   (pointer to lookup table of tangent values for CORDIC)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    l32r a13,.SAT_MAX_16  // Load maximum int16_t value
    l32r a14,.SAT_MIN_16  // Load minimum int16_t value

    l32r a11,.PLUS_PI     // Load +PI constant
    blt a2, a11,.sincos_n1  // Branch if angle < +PI
    add.n a2,a2,a14 // angle > +PI: angle = angle - 2*PI (a14 = -32768, so effectively subtracts -32768 which is add 32768, but wraps around for 2's complement overflow)
    j .sincos_q34         // Jump to quadrant determination
.sincos_n1:
    l32r a11,.MINUS_PI    // Load -PI constant
    bge a2, a11,.sincos_n2 // Branch if angle >= -PI
    add.n a2,a2,a13 // angle < -PI: angle = angle + 2*PI (a13 = 32767, similar wrapping logic applies)
.sincos_n2:

    l32r a11,.PLUS_PI2    // Load +PI/2 constant
    bge a2, a11, .sincos_q2 // Branch if angle >= PI/2 (Quadrant 2 or 3)
.sincos_q34:             // Angle is in Quadrant 3 or 4 (or Quadrant 1, handled here too)
    l32r a11,.MINUS_PI2   // Load -PI/2 constant
    blt a2, a11, .sincos_q3 // Branch if angle < -PI/2 (Quadrant 3)
    movi.n  a11, 0      // cor = 0 (no correction needed for Q1 or Q4) -> a11
    j .sincos_cordit11    // Jump to CORDIC calculation
.sincos_q2:             // Angle is in Quadrant 2
    l32r a11,.PLUS_PI   // cor = +PI -> a11 (correction factor)
    sub a2,a11,a2       // angle = PI - angle (fold angle into Q1)
    j .sincos_cordit11    // Jump to CORDIC calculation
.sincos_q3:             // Angle is in Quadrant 3
    l32r a11,.MINUS_PI  // cor = -PI -> a11 (correction factor)
    sub a2,a11,a2       // angle = -PI - angle (fold angle into Q1, result will be negative)
.sincos_cordit11:       // Start CORDIC calculation

    // sum -> a2 (input angle, potentially modified by quadrant folding), y -> a6 (sine output), x -> a7 (cosine output), xn -> a8 (temp x calc)
    movi.n  a6, 0       // y = 0 (initial sine value) -> a6
    l32r a7,.COS_FACTOR // x = initial gain factor for cosine -> a7

    // 0th CORDIC iteration
    l16si   a9, a5, 0   // Load tan[0] from lookup table -> a9
    bltz a2, .sincos_step0 // Branch if sum < 0
        // sum >= 0 (adjust rotation direction)
        sub a8,a7,a6    // xn = x - y (CORDIC x update)
        add.n a6,a6,a7  // y = y + x (CORDIC y update)
        sub a2,a2,a9  // sum -= tan[0] (adjust angle accumulator)
        j .sincos_step0_end // Jump to end of step 0
    .sincos_step0:
        // sum < 0 (adjust rotation direction)
        add.n a8,a7,a6  // xn = x + y (CORDIC x update)
        sub a6,a6,a7    // y = y - x (CORDIC y update)
        add.n a2,a2,a9  // sum += tan[0] (adjust angle accumulator)
    .sincos_step0_end:
        mov.n a7,a8     // x = xn (update x for next iteration)

    // 1st CORDIC iteration
    l16si a9, a5, 2     // Load tan[1] from lookup table -> a9
    srai a10, a7, 1      // a10 = x >> 1 (pre-shifted x for CORDIC)
    srai a12, a6, 1     // a12 = y >> 1 (pre-shifted y for CORDIC)
    bltz a2, .sincos_step1 // Branch if sum < 0
        // sum >= 0 (adjust rotation direction)
        sub a8,a7,a12    // xn = x - (y >> 1) (CORDIC x update)
        add.n a6,a6,a10  // y = y + (x >> 1) (CORDIC y update)
        sub a2,a2,a9  // sum -= tan[1] (adjust angle accumulator)
        j .sincos_step1_end // Jump to end of step 1
    .sincos_step1:
        // sum < 0 (adjust rotation direction)
        add.n a8,a7,a12  // xn = x + (y >> 1) (CORDIC x update)
        sub a6,a6,a10    // y = y - (x >> 1) (CORDIC y update)
        add.n a2,a2,a9  // sum += tan[1] (adjust angle accumulator)
    .sincos_step1_end:
        mov.n a7,a8    // x = xn (update x for next iteration)

    // 2nd CORDIC iteration
    l16si a9, a5, 4     // Load tan[2] from lookup table -> a9
    srai a10, a7, 2      // a10 = x >> 2 (pre-shifted x for CORDIC)
    srai a12, a6, 2     // a12 = y >> 2 (pre-shifted y for CORDIC)
    bltz a2, .sincos_step2 // Branch if sum < 0
        // sum >= 0 (adjust rotation direction)
        sub a8,a7,a12    // xn = x - (y >> 2) (CORDIC x update)
        add.n a6,a6,a10  // y = y + (x >> 2) (CORDIC y update)
        sub a2,a2,a9  // sum -= tan[2] (adjust angle accumulator)
        j .sincos_step2_end // Jump to end of step 2
    .sincos_step2:
        // sum < 0 (adjust rotation direction)
        add.n a8,a7,a12  // xn = x + (y >> 2) (CORDIC x update)
        sub a6,a6,a10    // y = y - (x >> 2) (CORDIC y update)
        add.n a2,a2,a9  // sum += tan[2] (adjust angle accumulator)
    .sincos_step2_end:
        mov.n a7,a8    // x = xn (update x for next iteration)

    // 3rd CORDIC iteration
    l16si a9, a5, 6     // Load tan[3] from lookup table -> a9
    srai a10, a7, 3      // a10 = x >> 3 (pre-shifted x for CORDIC)
    srai a12, a6, 3     // a12 = y >> 3 (pre-shifted y for CORDIC)
    bltz a2, .sincos_step3 // Branch if sum < 0
        // sum >= 0 (adjust rotation direction)
        sub a8,a7,a12    // xn = x - (y >> 3) (CORDIC x update)
        add.n a6,a6,a10  // y = y + (x >> 3) (CORDIC y update)
        sub a2,a2,a9  // sum -= tan[3] (adjust angle accumulator)
        j .sincos_step3_end // Jump to end of step 3
    .sincos_step3:
        // sum < 0 (adjust rotation direction)
        add.n a8,a7,a12  // xn = x + (y >> 3) (CORDIC x update)
        sub a6,a6,a10    // y = y - (x >> 3) (CORDIC y update)
        add.n a2,a2,a9  // sum += tan[3] (adjust angle accumulator)
    .sincos_step3_end:
        mov.n a7,a8    // x = xn (update x for next iteration)

    // 4th CORDIC iteration
    l16si a9, a5, 8     // Load tan[4] from lookup table -> a9
    srai a10, a7, 4      // a10 = x >> 4 (pre-shifted x for CORDIC)
    srai a12, a6, 4     // a12 = y >> 4 (pre-shifted y for CORDIC)
    bltz a2, .sincos_step4 // Branch if sum < 0
        // sum >= 0 (adjust rotation direction)
        sub a8,a7,a12    // xn = x - (y >> 4) (CORDIC x update)
        add.n a6,a6,a10  // y = y + (x >> 4) (CORDIC y update)
        sub a2,a2,a9  // sum -= tan[4] (adjust angle accumulator)
        j .sincos_step4_end // Jump to end of step 4
    .sincos_step4:
        // sum < 0 (adjust rotation direction)
        add.n a8,a7,a12  // xn = x + (y >> 4) (CORDIC x update)
        sub a6,a6,a10    // y = y - (x >> 4) (CORDIC y update)
        add.n a2,a2,a9  // sum += tan[4] (adjust angle accumulator)
    .sincos_step4_end:
        mov.n a7,a8    // x = xn (update x for next iteration)

    // 5th CORDIC iteration
    l16si a9, a5, 10     // Load tan[5] from lookup table -> a9
    srai a10, a7, 5      // a10 = x >> 5 (pre-shifted x for CORDIC)
    srai a12, a6, 5     // a12 = y >> 5 (pre-shifted y for CORDIC)
    bltz a2, .sincos_step5 // Branch if sum < 0
        // sum >= 0 (adjust rotation direction)
        sub a8,a7,a12    // xn = x - (y >> 5) (CORDIC x update)
        add.n a6,a6,a10  // y = y + (x >> 5) (CORDIC y update)
        sub a2,a2,a9  // sum -= tan[5] (adjust angle accumulator)
        j .sincos_step5_end // Jump to end of step 5
    .sincos_step5:
        // sum < 0 (adjust rotation direction)
        add.n a8,a7,a12  // xn = x + (y >> 5) (CORDIC x update)
        sub a6,a6,a10    // y = y - (x >> 5) (CORDIC y update)
        add.n a2,a2,a9  // sum += tan[5] (adjust angle accumulator)
    .sincos_step5_end:
        mov.n a7,a8    // x = xn (update x for next iteration)

    // 6th CORDIC iteration
    l16si a9, a5, 12     // Load tan[6] from lookup table -> a9
    srai a10, a7, 6      // a10 = x >> 6 (pre-shifted x for CORDIC)
    srai a12, a6, 6     // a12 = y >> 6 (pre-shifted y for CORDIC)
    bltz a2, .sincos_step6 // Branch if sum < 0
        // sum >= 0 (adjust rotation direction)
        sub a8,a7,a12    // xn = x - (y >> 6) (CORDIC x update)
        add.n a6,a6,a10  // y = y + (x >> 6) (CORDIC y update)
        sub a2,a2,a9  // sum -= tan[6] (adjust angle accumulator)
        j .sincos_step6_end // Jump to end of step 6
    .sincos_step6:
        // sum < 0 (adjust rotation direction)
        add.n a8,a7,a12  // xn = x + (y >> 6) (CORDIC x update)
        sub a6,a6,a10    // y = y - (x >> 6) (CORDIC y update)
        add.n a2,a2,a9  // sum += tan[6] (adjust angle accumulator)
    .sincos_step6_end:
        mov.n a7,a8    // x = xn (update x for next iteration)

    // 7th CORDIC iteration
    l16si a9, a5, 14     // Load tan[7] from lookup table -> a9
    srai a10, a7, 7      // a10 = x >> 7 (pre-shifted x for CORDIC)
    srai a12, a6, 7     // a12 = y >> 7 (pre-shifted y for CORDIC)
    bltz a2, .sincos_step7 // Branch if sum < 0
        // sum >= 0 (adjust rotation direction)
        sub a8,a7,a12    // xn = x - (y >> 7) (CORDIC x update)
        add.n a6,a6,a10  // y = y + (x >> 7) (CORDIC y update)
        sub a2,a2,a9  // sum -= tan[7] (adjust angle accumulator)
        j .sincos_step7_end // Jump to end of step 7
    .sincos_step7:
        // sum < 0 (adjust rotation direction)
        add.n a8,a7,a12  // xn = x + (y >> 7) (CORDIC x update)
        sub a6,a6,a10    // y = y - (x >> 7) (CORDIC y update)
        add.n a2,a2,a9  // sum += tan[7] (adjust angle accumulator)
    .sincos_step7_end:
        mov.n a7,a8    // x = xn (update x for next iteration)

    // 8th CORDIC iteration
    l16si a9, a5, 16     // Load tan[8] from lookup table -> a9
    srai a10, a7, 8      // a10 = x >> 8 (pre-shifted x for CORDIC)
    srai a12, a6, 8     // a12 = y >> 8 (pre-shifted y for CORDIC)
    bltz a2, .sincos_step8 // Branch if sum < 0
        // sum >= 0 (adjust rotation direction)
        sub a8,a7,a12    // xn = x - (y >> 8) (CORDIC x update)
        add.n a6,a6,a10  // y = y + (x >> 8) (CORDIC y update)
        sub a2,a2,a9  // sum -= tan[8] (adjust angle accumulator)
        j .sincos_step8_end // Jump to end of step 8
    .sincos_step8:
        // sum < 0 (adjust rotation direction)
        add.n a8,a7,a12  // xn = x + (y >> 8) (CORDIC x update)
        sub a6,a6,a10    // y = y - (x >> 8) (CORDIC y update)
        add.n a2,a2,a9  // sum += tan[8] (adjust angle accumulator)
    .sincos_step8_end:
        mov.n a7,a8    // x = xn (update x for next iteration)

    // 9th CORDIC iteration
    l16si a9, a5, 18     // Load tan[9] from lookup table -> a9
    srai a10, a7, 9      // a10 = x >> 9 (pre-shifted x for CORDIC)
    srai a12, a6, 9     // a12 = y >> 9 (pre-shifted y for CORDIC)
    bltz a2, .sincos_step9 // Branch if sum < 0
        // sum >= 0 (adjust rotation direction)
        sub a8,a7,a12    // xn = x - (y >> 9) (CORDIC x update)
        add.n a6,a6,a10  // y = y + (x >> 9) (CORDIC y update)
        sub a2,a2,a9  // sum -= tan[9] (adjust angle accumulator)
        j .sincos_step9_end // Jump to end of step 9
    .sincos_step9:
        // sum < 0 (adjust rotation direction)
        add.n a8,a7,a12  // xn = x + (y >> 9) (CORDIC x update)
        sub a6,a6,a10    // y = y - (x >> 9) (CORDIC y update)
        add.n a2,a2,a9  // sum += tan[9] (adjust angle accumulator)
    .sincos_step9_end:
        mov.n a7,a8    // x = xn (update x for next iteration)

    // 10th CORDIC iteration
    l16si a9, a5, 20     // Load tan[10] from lookup table -> a9
    srai a10, a7, 10      // a10 = x >> 10 (pre-shifted x for CORDIC)
    srai a12, a6, 10     // a12 = y >> 10 (pre-shifted y for CORDIC)
    bltz a2, .sincos_step10 // Branch if sum < 0
        // sum >= 0 (adjust rotation direction)
        sub a8,a7,a12    // xn = x - (y >> 10) (CORDIC x update)
        add.n a6,a6,a10  // y = y + (x >> 10) (CORDIC y update)
        sub a2,a2,a9  // sum -= tan[10] (adjust angle accumulator)
        j .sincos_step10_end // Jump to end of step 10
    .sincos_step10:
        // sum < 0 (adjust rotation direction)
        add.n a8,a7,a12  // xn = x + (y >> 10) (CORDIC x update)
        sub a6,a6,a10    // y = y - (x >> 10) (CORDIC y update)
        add.n a2,a2,a9  // sum += tan[10] (adjust angle accumulator)
    .sincos_step10_end:
        mov.n a7,a8    // x = xn (update x for next iteration)

    // 11th CORDIC iteration
    l16si a9, a5, 22     // Load tan[11] from lookup table -> a9
    srai a10, a7, 11      // a10 = x >> 11 (pre-shifted x for CORDIC)
    srai a12, a6, 11     // a12 = y >> 11 (pre-shifted y for CORDIC)
    bltz a2, .sincos_step11 // Branch if sum < 0
        // sum >= 0 (adjust rotation direction)
        sub a8,a7,a12    // xn = x - (y >> 11) (CORDIC x update)
        add.n a6,a6,a10  // y = y + (x >> 11) (CORDIC y update)
        sub a2,a2,a9  // sum -= tan[11] (adjust angle accumulator)
        j .sincos_step11_end // Jump to end of step 11
    .sincos_step11:
        // sum < 0 (adjust rotation direction)
        add.n a8,a7,a12  // xn = x + (y >> 11) (CORDIC x update)
        sub a6,a6,a10    // y = y - (x >> 11) (CORDIC y update)
        add.n a2,a2,a9  // sum += tan[11] (adjust angle accumulator)
    .sincos_step11_end:
        mov.n a7,a8    // x = xn (update x for next iteration)

    // 12th CORDIC iteration
    l16si a9, a5, 24     // Load tan[12] from lookup table -> a9
    srai a10, a7, 12      // a10 = x >> 12 (pre-shifted x for CORDIC)
    srai a12, a6, 12     // a12 = y >> 12 (pre-shifted y for CORDIC)
    bltz a2, .sincos_step12 // Branch if sum < 0
        // sum >= 0 (adjust rotation direction)
        sub a8,a7,a12    // xn = x - (y >> 12) (CORDIC x update)
        add.n a6,a6,a10  // y = y + (x >> 12) (CORDIC y update)
        sub a2,a2,a9  // sum -= tan[12] (adjust angle accumulator)
        j .sincos_step12_end // Jump to end of step 12
    .sincos_step12:
        // sum < 0 (adjust rotation direction)
        add.n a8,a7,a12  // xn = x + (y >> 12) (CORDIC x update)
        sub a6,a6,a10    // y = y - (x >> 12) (CORDIC y update)
        add.n a2,a2,a9  // sum += tan[12] (adjust angle accumulator)
    .sincos_step12_end:
        mov.n a7,a8    // x = xn (update x for next iteration)

    // 13th CORDIC iteration (last one, no lookup table access)
    srai a10, a7, 13      // a10 = x >> 13 (pre-shifted x for CORDIC)
    srai a12, a6, 13     // a12 = y >> 13 (pre-shifted y for CORDIC)
    bltz a2, .sincos_step13 // Branch if sum < 0
        // sum >= 0 (adjust rotation direction)
        sub a8,a7,a12    // xn = x - (y >> 13) (CORDIC x update)
        add.n a6,a6,a10  // y = y + (x >> 13) (CORDIC y update)
        j .sincos_step13_end // Jump to end of step 13
    .sincos_step13:
        // sum < 0 (adjust rotation direction)
        add.n a8,a7,a12  // xn = x + (y >> 13) (CORDIC x update)
        sub a6,a6,a10    // y = y - (x >> 13) (CORDIC y update)
    .sincos_step13_end:
        // x = xn -> a8 (final x value is in a8)

    beqz a11,.sincos_end  // Branch if cor (a11) is 0 (no quadrant correction needed)
    neg a8,a8             // Negate x if correction was applied (for Q2/Q3)
.sincos_end:

    // Saturate results to fit in int16_t range
    //l32r a13,.SAT_MAX_16 (a13 already loaded)
    min a8,a8,a13  // Clamp x (cos) to max int16_t
    min a6,a6,a13  // Clamp y (sin) to max int16_t
    //l32r a14,.SAT_MIN_16 (a14 already loaded)
    max a8,a8,a14  // Clamp x (cos) to min int16_t
    max a6,a6,a14  // Clamp y (sin) to min int16_t
    
    s16i a8,a3,0 // Store saturated cosine value to address pointed by a3
    s16i a6,a4,0 // Store saturated sine value to address pointed by a4
    retw.n       // Return from function (windowed return)

///////////////////////////////////////////////////////////////////////////////////////////////////////
