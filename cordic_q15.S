    .text
    .align  4
	.literal_position
	.literal	.PLUS_PI, 16383
	.literal	.MINUS_PI, -16383
	.literal	.PLUS_PI2, 8191
	.literal	.MINUS_PI2, -8191
//    .literal	.COS_FACTOR, 20860
    .literal	.COS_FACTOR, 19894
	.literal	.SAT_MAX_16, 32767
	.literal	.SAT_MIN_16, -32768
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   atan2_q15_s3,@function
    .align  4
    .global atan2_q15_s3
atan2_q15_s3: 
// y   - a2
// x   - a3
// &tan - a4
    entry	a1, 16

    abs 	a5, a3       // x1 = abs(x)	
    movi.n  a6, 0        // sum = 0

    // 0 step
    l16si   a8, a4, 0   // a8 = tan[0]
    bgez a2, .atan2_step0
        // y < 0
        sub a7,a5,a2    // xn = x1 - y
        add.n a2,a2,a5  // y = y + x1
        sub a6,a6,a8  // sum -= tan[0]
        j .atan2_step0_end
    .atan2_step0:
        // y >= 0
        add.n a7,a5,a2  // xn = x1 + y
        sub a2,a2,a5    // y = y - x1
        add.n a6,a6,a8  // sum += tan[0]
    .atan2_step0_end:
        beqz a2,.atan2_step12_end
        mov.n a5,a7     //  x1 = xn

    // 1 step
    l16si a8, a4, 2     // a8 = tan[1]
    srai a9, a5, 1      // a9 = (x1 >> 1)    
    srai a10, a2, 1     // a10 = (y >> 1)    
    bgez a2, .atan2_step1
        // y < 0
        sub a7,a5,a10    // xn = x1 - (y >> 1)
        add.n a2,a2,a9  // y = y + (x1 >> 1) 
        sub a6,a6,a8  // sum -= tan[1]
        j .atan2_step1_end
    .atan2_step1:
        // y >= 0
        add.n a7,a5,a10  // xn = x1 + (y >> 1)
        sub a2,a2,a9    // y = y - (x1 >> 1) 
        add.n a6,a6,a8  // sum += tan[1]
    .atan2_step1_end:
        beqz a2,.atan2_step12_end
        mov.n a5,a7    //  x1 = xn

    // 2 step
    l16si a8, a4, 4     // a8 = tan[2]
    srai a9, a5, 2      // a9 = (x1 >> 2)    
    srai a10, a2, 2     // a10 = (y >> 2)    
    bgez a2, .atan2_step2
        // y < 0
        sub a7,a5,a10    // xn = x1 - (y >> 2)
        add.n a2,a2,a9  // y = y + (x1 >> 2) 
        sub a6,a6,a8  // sum -= tan[2]
        j .atan2_step2_end
    .atan2_step2:
        // y >= 0
        add.n a7,a5,a10  // xn = x1 + (y >> 2)
        sub a2,a2,a9    // y = y - (x1 >> 2) 
        add.n a6,a6,a8  // sum += tan[2]
    .atan2_step2_end:
        beqz a2,.atan2_step12_end
        mov.n a5,a7    //  x1 = xn

    // 3 step
    l16si a8, a4, 6     // a8 = tan[3]
    srai a9, a5, 3      // a9 = (x1 >> 3)    
    srai a10, a2, 3     // a10 = (y >> 3)    
    bgez a2, .atan2_step3
        // y < 0
        sub a7,a5,a10    // xn = x1 - (y >> 3)
        add.n a2,a2,a9  // y = y + (x1 >> 3) 
        sub a6,a6,a8  // sum -= tan[3]
        j .atan2_step3_end
    .atan2_step3:
        // y >= 0
        add.n a7,a5,a10  // xn = x1 + (y >> 3)
        sub a2,a2,a9    // y = y - (x1 >> 3) 
        add.n a6,a6,a8  // sum += tan[3]
    .atan2_step3_end:
        beqz a2,.atan2_step12_end
        mov.n a5,a7    //  x1 = xn

    // 4 step
    l16si a8, a4, 8     // a8 = tan[4]
    srai a9, a5, 4      // a9 = (x1 >> 4)    
    srai a10, a2, 4     // a10 = (y >> 4)    
    bgez a2, .atan2_step4
        // y < 0
        sub a7,a5,a10    // xn = x1 - (y >> 4)
        add.n a2,a2,a9  // y = y + (x1 >> 4) 
        sub a6,a6,a8  // sum -= tan[4]
        j .atan2_step4_end
    .atan2_step4:
        // y >= 0
        add.n a7,a5,a10  // xn = x1 + (y >> 4)
        sub a2,a2,a9    // y = y - (x1 >> 4) 
        add.n a6,a6,a8  // sum += tan[4]
    .atan2_step4_end:
        beqz a2,.atan2_step12_end
        mov.n a5,a7    //  x1 = xn

    // 5 step
    l16si a8, a4, 10     // a8 = tan[5]
    srai a9, a5, 5      // a9 = (x1 >> 5)    
    srai a10, a2, 5     // a10 = (y >> 5)    
    bgez a2, .atan2_step5
        // y < 0
        sub a7,a5,a10    // xn = x1 - (y >> 5)
        add.n a2,a2,a9  // y = y + (x1 >> 5) 
        sub a6,a6,a8  // sum -= tan[5]
        j .atan2_step5_end
    .atan2_step5:
        // y >= 0
        add.n a7,a5,a10  // xn = x1 + (y >> 5)
        sub a2,a2,a9    // y = y - (x1 >> 5) 
        add.n a6,a6,a8  // sum += tan[5]
    .atan2_step5_end:
        beqz a2,.atan2_step12_end
        mov.n a5,a7    //  x1 = xn

    // 6 step
    l16si a8, a4, 12     // a8 = tan[6]
    srai a9, a5, 6      // a9 = (x1 >> 6)    
    srai a10, a2, 6     // a10 = (y >> 6)    
    bgez a2, .atan2_step6
        // y < 0
        sub a7,a5,a10    // xn = x1 - (y >> 6)
        add.n a2,a2,a9  // y = y + (x1 >>6) 
        sub a6,a6,a8  // sum -= tan[6]
        j .atan2_step6_end
    .atan2_step6:
        // y >= 0
        add.n a7,a5,a10  // xn = x1 + (y >> 6)
        sub a2,a2,a9    // y = y - (x1 >> 6) 
        add.n a6,a6,a8  // sum += tan[6]
    .atan2_step6_end:
        beqz a2,.atan2_step12_end
        mov.n a5,a7    //  x1 = xn

    // 7 step
    l16si a8, a4, 14     // a8 = tan[7]
    srai a9, a5, 7      // a9 = (x1 >> 7)    
    srai a10, a2, 7     // a10 = (y >> 7)    
    bgez a2, .atan2_step7
        // y < 0
        sub a7,a5,a10    // xn = x1 - (y >> 7)
        add.n a2,a2,a9  // y = y + (x1 >>7) 
        sub a6,a6,a8  // sum -= tan[7]
        j .atan2_step7_end
    .atan2_step7:
        // y >= 0
        add.n a7,a5,a10  // xn = x1 + (y >> 7)
        sub a2,a2,a9    // y = y - (x1 >> 7) 
        add.n a6,a6,a8  // sum += tan[7]
    .atan2_step7_end:
        beqz a2,.atan2_step12_end
        mov.n a5,a7    //  x1 = xn

    // 8 step
    l16si a8, a4, 16     // a8 = tan[8]
    srai a9, a5, 8      // a9 = (x1 >> 8)    
    srai a10, a2, 8     // a10 = (y >> 8)    
    bgez a2, .atan2_step8
        // y < 0
        sub a7,a5,a10    // xn = x1 - (y >> 8)
        add.n a2,a2,a9  // y = y + (x1 >>8) 
        sub a6,a6,a8  // sum -= tan[8]
        j .atan2_step8_end
    .atan2_step8:
        // y >= 0
        add.n a7,a5,a10  // xn = x1 + (y >> 8)
        sub a2,a2,a9    // y = y - (x1 >> 8) 
        add.n a6,a6,a8  // sum += tan[8]
    .atan2_step8_end:
        beqz a2,.atan2_step12_end
        mov.n a5,a7    //  x1 = xn

    // 9 step
    l16si a8, a4, 18     // a8 = tan[9]
    srai a9, a5, 9      // a9 = (x1 >> 9)    
    srai a10, a2, 9     // a10 = (y >> 9)    
    bgez a2, .atan2_step9
        // y < 0
        sub a7,a5,a10    // xn = x1 - (y >> 9)
        add.n a2,a2,a9  // y = y + (x1 >>9) 
        sub a6,a6,a8  // sum -= tan[9]
        j .atan2_step9_end
    .atan2_step9:
        // y >= 0
        add.n a7,a5,a10  // xn = x1 + (y >> 9)
        sub a2,a2,a9    // y = y - (x1 >> 9) 
        add.n a6,a6,a8  // sum += tan[9]
    .atan2_step9_end:
        beqz a2,.atan2_step12_end
        mov.n a5,a7    //  x1 = xn

    // 10 step
    l16si a8, a4, 20     // a8 = tan[10]
    srai a9, a5, 10      // a9 = (x1 >> 10)    
    srai a10, a2, 10     // a10 = (y >> 10)    
    bgez a2, .atan2_step10
        // y < 0
        sub a7,a5,a10    // xn = x1 - (y >> 10)
        add.n a2,a2,a9  // y = y + (x1 >>10) 
        sub a6,a6,a8  // sum -= tan[10]
        j .atan2_step10_end
    .atan2_step10:
        // y >= 0
        add.n a7,a5,a10  // xn = x1 + (y >> 10)
        sub a2,a2,a9    // y = y - (x1 >> 10) 
        add.n a6,a6,a8  // sum += tan[10]
    .atan2_step10_end:
        beqz a2,.atan2_step12_end
        mov.n a5,a7    //  x1 = xn

    // 11 step
    l16si a8, a4, 22     // a8 = tan[11]
    srai a9, a5, 11      // a9 = (x1 >> 11)    
    srai a10, a2, 11     // a10 = (y >> 11)    
    bgez a2, .atan2_step11
        // y < 0
        sub a7,a5,a10    // xn = x1 - (y >> 11)
        add.n a2,a2,a9  // y = y + (x1 >>11) 
        sub a6,a6,a8  // sum -= tan[11]
        j .atan2_step11_end
    .atan2_step11:
        // y >= 0
        add.n a7,a5,a10  // xn = x1 + (y >> 11)
        sub a2,a2,a9    // y = y - (x1 >> 11) 
        add.n a6,a6,a8  // sum += tan[11]
    .atan2_step11_end:
        beqz a2,.atan2_step12_end

    // 12 step
    l16si a8, a4, 24     // a8 = tan[12]
    bgez a2, .atan2_step12
        // y < 0
        sub a6,a6,a8  // sum -= tan[12]
        j .atan2_step12_end
    .atan2_step12:
        // y >= 0
        add.n a6,a6,a8  // sum += tan[12]
    .atan2_step12_end:
   
    bgez a3,.end_atan2_q15_s3
        // x < 0
        bgez a6,.atan2_last_check
            // sum < 0
            l32r a7,.MINUS_PI
            j .end1_atan2_q15_s3
        .atan2_last_check:
            // sum >= 0
            l32r a7,.PLUS_PI
        .end1_atan2_q15_s3:
        sub a2,a7,a6
        retw.n
    .end_atan2_q15_s3:
    mov.n a2,a6
    retw.n
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   sincos_q15_s3,@function
    .align  4
    .global sincos_q15_s3
sincos_q15_s3: 
// angle   - a2
// &cos   - a3
// &sin   - a4
// &tan - a5
    entry	a1, 16

    l32r a13,.SAT_MAX_16
    l32r a14,.SAT_MIN_16

    l32r a11,.PLUS_PI
    blt a2, a11,.sincos_n1
    add.n a2,a2,a14 // >Pi
    j .sincos_q34
.sincos_n1:
    l32r a11,.MINUS_PI
    bge a2, a11,.sincos_n2
    add.n a2,a2,a13 // <-Pi
.sincos_n2:

    l32r a11,.PLUS_PI2
    bge a2, a11, .sincos_q2
.sincos_q34:
    l32r a11,.MINUS_PI2
    blt a2, a11, .sincos_q3
    movi.n  a11, 0      // cor = 0              -> a11
    j .sincos_cordit11
.sincos_q2:             // if angle >= pi/2
    l32r a11,.PLUS_PI   // cor = 18383          -> a11
    sub a2,a11,a2       // angle = cor - angle  -> a2
    j .sincos_cordit11
.sincos_q3:             // if angle <= pi/2
    l32r a11,.MINUS_PI  // cor = -18383         -> a11
    sub a2,a11,a2       // angle = cor - angle  -> a2
.sincos_cordit11:

    // sum -> a2, y -> a6, x -> a7, xn -> a8
    movi.n  a6, 0       // y = 0                -> a6
    l32r a7,.COS_FACTOR // x = gain1            -> a7

    // 0 step
    l16si   a9, a5, 0   // tan[0]               -> a9
    bltz a2, .sincos_step0
        // sum > 0
        sub a8,a7,a6    // xn = x - y
        add.n a6,a6,a7  // y = y + x
        sub a2,a2,a9  // sum -= tan[0]
        j .sincos_step0_end
    .sincos_step0:
        // sum <= 0
        add.n a8,a7,a6  // xn = x + y
        sub a6,a6,a7    // y = y - x
        add.n a2,a2,a9  // sum += tan[0]
    .sincos_step0_end:
        mov.n a7,a8     //  x = xn

    // 1 step
    l16si a9, a5, 2     // tan[1]               -> a9
    srai a10, a7, 1      // a10 = (x >> 1)    
    srai a12, a6, 1     // a12 = (y >> 1)    
    bltz a2, .sincos_step1
        // sum >= 0
        sub a8,a7,a12    // xn = x - (y >> 1)
        add.n a6,a6,a10  // y = y + (x >> 1) 
        sub a2,a2,a9  // sum -= tan[1]
        j .sincos_step1_end
    .sincos_step1:
        // sum < 0
        add.n a8,a7,a12  // xn = x + (y >> 1)
        sub a6,a6,a10    // y = y - (x >> 1) 
        add.n a2,a2,a9  // sum += tan[1]
    .sincos_step1_end:
        mov.n a7,a8    //  x = xn

    // 2 step
    l16si a9, a5, 4     // tan[2]               -> a9
    srai a10, a7, 2      // a10 = (x >> 2)    
    srai a12, a6, 2     // a12 = (y >> 2)    
    bltz a2, .sincos_step2
        // sum >= 0
        sub a8,a7,a12    // xn = x - (y >> 2)
        add.n a6,a6,a10  // y = y + (x >> 2) 
        sub a2,a2,a9  // sum -= tan[]
        j .sincos_step2_end
    .sincos_step2:
        // sum < 0
        add.n a8,a7,a12  // xn = x + (y >> 2)
        sub a6,a6,a10    // y = y - (x >> 2) 
        add.n a2,a2,a9  // sum += tan[]
    .sincos_step2_end:
        mov.n a7,a8    //  x = xn

    // 3 step
    l16si a9, a5, 6     // tan[3]               -> a9
    srai a10, a7, 3      // a10 = (x >> 3)    
    srai a12, a6, 3     // a12 = (y >> 3)    
    bltz a2, .sincos_step3
        // sum >= 0
        sub a8,a7,a12    // xn = x - (y >> 3)
        add.n a6,a6,a10  // y = y + (x >> 3) 
        sub a2,a2,a9  // sum -= tan[]
        j .sincos_step3_end
    .sincos_step3:
        // sum < 0
        add.n a8,a7,a12  // xn = x + (y >> 3)
        sub a6,a6,a10    // y = y - (x >> 3) 
        add.n a2,a2,a9  // sum += tan[]
    .sincos_step3_end:
        mov.n a7,a8    //  x = xn

    // 4 step
    l16si a9, a5, 8     // tan[4]               -> a9
    srai a10, a7, 4      // a10 = (x >> 4)    
    srai a12, a6, 4     // a12 = (y >> 4)    
    bltz a2, .sincos_step4
        // sum >= 0
        sub a8,a7,a12    // xn = x - (y >> 4)
        add.n a6,a6,a10  // y = y + (x >> 4) 
        sub a2,a2,a9  // sum -= tan[]
        j .sincos_step4_end
    .sincos_step4:
        // sum < 0
        add.n a8,a7,a12  // xn = x + (y >> 4)
        sub a6,a6,a10    // y = y - (x >> 4) 
        add.n a2,a2,a9  // sum += tan[]
    .sincos_step4_end:
        mov.n a7,a8    //  x = xn

    // 5 step
    l16si a9, a5, 10     // tan[5]               -> a9
    srai a10, a7, 5      // a10 = (x >> 5)    
    srai a12, a6, 5     // a12 = (y >> 5)    
    bltz a2, .sincos_step5
        // sum >= 0
        sub a8,a7,a12    // xn = x - (y >> 5)
        add.n a6,a6,a10  // y = y + (x >> 5) 
        sub a2,a2,a9  // sum -= tan[]
        j .sincos_step5_end
    .sincos_step5:
        // sum < 0
        add.n a8,a7,a12  // xn = x + (y >> 5)
        sub a6,a6,a10    // y = y - (x >> 5) 
        add.n a2,a2,a9  // sum += tan[]
    .sincos_step5_end:
        mov.n a7,a8    //  x = xn

    // 6 step
    l16si a9, a5, 12     // tan[6]               -> a9
    srai a10, a7, 6      // a10 = (x >> 6)    
    srai a12, a6, 6     // a12 = (y >> 6)    
    bltz a2, .sincos_step6
        // sum >= 0
        sub a8,a7,a12    // xn = x - (y >> 6)
        add.n a6,a6,a10  // y = y + (x >> 6) 
        sub a2,a2,a9  // sum -= tan[]
        j .sincos_step6_end
    .sincos_step6:
        // sum < 0
        add.n a8,a7,a12  // xn = x + (y >> 6)
        sub a6,a6,a10    // y = y - (x >> 6) 
        add.n a2,a2,a9  // sum += tan[]
    .sincos_step6_end:
        mov.n a7,a8    //  x = xn

    // 7 step
    l16si a9, a5, 14     // tan[7]               -> a9
    srai a10, a7, 7      // a10 = (x >> 7)    
    srai a12, a6, 7     // a12 = (y >> 7)    
    bltz a2, .sincos_step7
        // sum >= 0
        sub a8,a7,a12    // xn = x - (y >> 7)
        add.n a6,a6,a10  // y = y + (x >> 7) 
        sub a2,a2,a9  // sum -= tan[]
        j .sincos_step7_end
    .sincos_step7:
        // sum < 0
        add.n a8,a7,a12  // xn = x + (y >> 7)
        sub a6,a6,a10    // y = y - (x >> 7) 
        add.n a2,a2,a9  // sum += tan[]
    .sincos_step7_end:
        mov.n a7,a8    //  x = xn

    // 8 step
    l16si a9, a5, 16     // tan[8]               -> a9
    srai a10, a7, 8      // a10 = (x >> 8)    
    srai a12, a6, 8     // a12 = (y >> 8)    
    bltz a2, .sincos_step8
        // sum >= 0
        sub a8,a7,a12    // xn = x - (y >> 8)
        add.n a6,a6,a10  // y = y + (x >> 8) 
        sub a2,a2,a9  // sum -= tan[]
        j .sincos_step8_end
    .sincos_step8:
        // sum < 0
        add.n a8,a7,a12  // xn = x + (y >> 8)
        sub a6,a6,a10    // y = y - (x >> 8) 
        add.n a2,a2,a9  // sum += tan[]
    .sincos_step8_end:
        mov.n a7,a8    //  x = xn

    // 9 step
    l16si a9, a5, 18     // tan[9]               -> a9
    srai a10, a7, 9      // a10 = (x >> 9)    
    srai a12, a6, 9     // a12 = (y >> 9)    
    bltz a2, .sincos_step9
        // sum >= 0
        sub a8,a7,a12    // xn = x - (y >> 9)
        add.n a6,a6,a10  // y = y + (x >> 9) 
        sub a2,a2,a9  // sum -= tan[]
        j .sincos_step9_end
    .sincos_step9:
        // sum < 0
        add.n a8,a7,a12  // xn = x + (y >> 9)
        sub a6,a6,a10    // y = y - (x >> 9) 
        add.n a2,a2,a9  // sum += tan[]
    .sincos_step9_end:
        mov.n a7,a8    //  x = xn

    // 10 step
    l16si a9, a5, 20     // tan[10]               -> a9
    srai a10, a7, 10      // a10 = (x >> 10)    
    srai a12, a6, 10     // a12 = (y >> 10)    
    bltz a2, .sincos_step10
        // sum >= 0
        sub a8,a7,a12    // xn = x - (y >> 10)
        add.n a6,a6,a10  // y = y + (x >> 10) 
        sub a2,a2,a9  // sum -= tan[]
        j .sincos_step10_end
    .sincos_step10:
        // sum < 0
        add.n a8,a7,a12  // xn = x + (y >> 10)
        sub a6,a6,a10    // y = y - (x >> 10) 
        add.n a2,a2,a9  // sum += tan[]
    .sincos_step10_end:
        mov.n a7,a8    //  x = xn

    // 11 step
    l16si a9, a5, 22     // tan[11]               -> a9
    srai a10, a7, 11      // a10 = (x >> 11)    
    srai a12, a6, 11     // a12 = (y >> 11)    
    bltz a2, .sincos_step11
        // sum >= 0
        sub a8,a7,a12    // xn = x - (y >> 11)
        add.n a6,a6,a10  // y = y + (x >> 11) 
        sub a2,a2,a9  // sum -= tan[]
        j .sincos_step11_end
    .sincos_step11:
        // sum < 0
        add.n a8,a7,a12  // xn = x + (y >> 11)
        sub a6,a6,a10    // y = y - (x >> 11) 
        add.n a2,a2,a9  // sum += tan[]
    .sincos_step11_end:
        mov.n a7,a8    //  x = xn

    // 12 step
    l16si a9, a5, 24     // tan[12]               -> a9
    srai a10, a7, 12      // a10 = (x >> 12)    
    srai a12, a6, 12     // a12 = (y >> 12)    
    bltz a2, .sincos_step12
        // sum >= 0
        sub a8,a7,a12    // xn = x - (y >> 12)
        add.n a6,a6,a10  // y = y + (x >> 12) 
        sub a2,a2,a9  // sum -= tan[]
        j .sincos_step12_end
    .sincos_step12:
        // sum < 0
        add.n a8,a7,a12  // xn = x + (y >> 12)
        sub a6,a6,a10    // y = y - (x >> 12) 
        add.n a2,a2,a9  // sum += tan[]
    .sincos_step12_end:
        mov.n a7,a8    //  x = xn

    // 13 step
    srai a10, a7, 13      // a10 = (x >> 13)    
    srai a12, a6, 13     // a12 = (y >> 13)    
    bltz a2, .sincos_step13
        // sum >= 0
        sub a8,a7,a12    // xn = x - (y >> 13)
        add.n a6,a6,a10  // y = y + (x >> 13) 
        j .sincos_step13_end
    .sincos_step13:
        // sum < 0
        add.n a8,a7,a12  // xn = x + (y >> 13)
        sub a6,a6,a10    // y = y - (x >> 13) 
    .sincos_step13_end:
        //  x = xn -> a8

    beqz a11,.sincos_end
    neg a8,a8
.sincos_end:

    // saturate
    //l32r a13,.SAT_MAX_16
    min a8,a8,a13
    min a6,a6,a13
    //l32r a14,.SAT_MIN_16
    max a8,a8,a14
    max a6,a6,a14
    
    s16i a8,a3,0 //cos
    s16i a6,a4,0 //sin
    retw.n
///////////////////////////////////////////////////////////////////////////////////////////////////////
