
    .text
    .align  4
	.literal_position
	.literal	.SAT_MAX_30, 1073692671
	.literal	.SAT_MIN_30, -1073692672
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   fir_16_16_q15_pie,@function
    .align  4
    .global fir_16_16_q15_pie
fir_16_16_q15_pie: 
// in   - a2
// k    - a3
// ksize- a4
// out  - a5
// size - a6
    entry	a1, 16

    srli 	a4, a4, 3           //  ksize/8 	
    slli 	a6, a6, 1           //  size*2 	
    add.n 	a6, a6, a5          //  a6 = out + size*2 end of out	

    EE.VLD.128.IP  q0, a2, 16   // q0 = in[0 .. 7]
    MOV.N a7, a2                // a7 = &in[8]
    MOV.N a8, a3                // a8 = k
    MOVI.N  a9, 15              // later shift right 15

.big_loop:
    EE.ZERO.QACC                        // QACC = 0
    loopnez a4, .loop_end_fir_pie
        EE.VLD.128.IP  q1, a2, 16       // q1 = in[i+8 .. i+15]
        EE.VLD.128.IP  q5, a3, 16       // q5 = k[i .. i+7]

        EE.VSMULAS.S16.QACC q0, q5, 0   // QACC[0..7] += in[i .. i+7] * k[i]
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+1 .. i+8], q1 = in[i+9 .. i+15 0]
        EE.VSMULAS.S16.QACC q0, q5, 1   // QACC[0..7] += in[i+1 .. i+8] * k[i+1]
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+2 .. i+9], q1 = in[i+10 .. i+15 0 0]
        EE.VSMULAS.S16.QACC q0, q5, 2   // QACC[0..7] += in[i+2 .. i+9] * k[i+2]
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+3 .. i+10], q1 = in[i+11 .. i+15, 0 .. 0]
        EE.VSMULAS.S16.QACC q0, q5, 3   // QACC[0..7] += in[i+3 .. i+10] * k[i+3]
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+4 .. i+11], q1 = in[i+12 .. i+15, 0 .. 0]
        EE.VSMULAS.S16.QACC q0, q5, 4   // QACC[0..7] += in[i+4 .. i+11] * k[i+4]
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+5 .. i+12], q1 = in[i+13 .. i+15, 0 .. 0]
        EE.VSMULAS.S16.QACC q0, q5, 5   // QACC[0..7] += in[i+5 .. i+12] * k[i+5]
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+6 .. i+13], q1 = in[i+14 i+15, 0 .. 0]
        EE.VSMULAS.S16.QACC q0, q5, 6   // QACC[0..7] += in[i+6 .. i+13] * k[i+6]
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+7 .. i+14], q1 = in[i+15, 0 .. 0]
        EE.VSMULAS.S16.QACC q0, q5, 7   // QACC[0..7] += in[i+7 .. i+14] * k[i+7]
        EE.SRCI.2Q q1, q0, 1            // q0 = in[i+8 .. i+15], q1 = 0
    .loop_end_fir_pie:                  // loop ksize/8
    EE.SRCMB.S16.QACC q1, a9, 0         // q1 = saturate(QACC[0..7] >> 15)
    EE.VST.128.IP q1, a5, 16            // out[i..i+7] = q1 
    BEQ a5,a6,.end_fir_16_16_q15_pie    // check end of out 

    EE.VLD.128.IP  q0, a7, 16       // q0 = in[i .. i+7]
    MOV.N a2, a7                    // a7 = &in[i+8]
    MOV.N a3, a8                    // a3 = k
    J .big_loop

.end_fir_16_16_q15_pie:
    retw.n
