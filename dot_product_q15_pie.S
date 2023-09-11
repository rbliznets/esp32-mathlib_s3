    .text
    .align  4
	.literal_position
	.literal	.SAT_MAX_30, 1073741823
	.literal	.SAT_MIN_30, -1073741824
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   dot_product_q15_16_16_pie,@function
    .align  4
    .global dot_product_q15_16_16_pie
dot_product_q15_16_16_pie: 
// in1   - a2
// in2   - a3
// size  - a4
    entry	a1, 16

    srli 	a4, a4, 3           //  size/8 	
    ee.zero.accx                //  accx = 0 	
    ee.vld.128.ip   q0, a2, 16  //  in1[0..7] -> q0	
    ee.vld.128.ip   q1, a3, 16  //  in2[0..7] -> q1	

    loopnez a4, .loop_end_dot_product_q15_16_16_pie
        ee.vmulas.s16.accx.ld.ip    q0, a2, 16, q0, q1  // accx += dotproduct(q0,q1);  in1[n..n+7] -> q0
        ee.vld.128.ip               q1, a3, 16          // in2[n..n+7] -> q1
    .loop_end_dot_product_q15_16_16_pie:

    movi.n	a5, 0
    EE.SRS.ACCX a2, a5, 0  // a2 = sat(accx)

    l32r	a3, .SAT_MAX_30
    min     a2, a2, a3          // if( a2 > 1073741823) a2 = 1073741823
    l32r	a4, .SAT_MIN_30
    max     a2, a2, a4          // if( a2 < -1073741824) a2 = -1073741824
    SRAI    a2, a2, 15          // a2 = a2 >> 15

    retw.n
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   dot_product_q15_1_16_pie,@function
    .align  4
    .global dot_product_q15_1_16_pie
dot_product_q15_1_16_pie: 
// in1   - a2
// in2   - a3
// size  - a4
    entry	a1, 16

    srli 	a4, a4, 3           //  size/8 	

    ee.zero.accx                //  accx = 0 	
    EE.LD.128.USAR.IP   q2, a2, 16 
    ee.vld.128.ip   q1, a3, 16  //  in2[0..7] -> q1	

    loopnez a4, .loop_end_dot_product_q15_1_16_pie
        EE.LD.128.USAR.IP   q3, a2, 16 
        EE.SRC.Q.QUP        q0, q2, q3  // in1[0..7] -> q0
        ee.vmulas.s16.accx.ld.ip    q1, a3, 16, q0, q1  // accx += dotproduct(q0,q1);  in2[n..n+7] -> q1
    .loop_end_dot_product_q15_1_16_pie:

    movi.n	a5, 0
    EE.SRS.ACCX a2, a5, 0  // a2 = sat(accx)

    l32r	a3, .SAT_MAX_30
    min     a2, a2, a3          // if( a2 > 1073741823) a2 = 1073741823
    l32r	a4, .SAT_MIN_30
    max     a2, a2, a4          // if( a2 < -1073741824) a2 = -1073741824
    SRAI    a2, a2, 15          // a2 = a2 >> 15

    retw.n
