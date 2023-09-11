    .type   arg_16_q15_pie,@function
    .align  4
    .global arg_16_q15_pie
arg_16_q15_pie: 
// in   - a2
// out  - a3
// size  - a4
// tan - a5
    entry	a1, 16
    
    MOVI.N a9, 13
    MOVI.N a7, 0
    MOVI.N a6, -1
    EE.ZERO.Q q7
    MOV.N a10, a5

    slli 	a4, a4, 1           //  size*2 	
    add.n 	a4, a4, a3          //  a4 = out + size*2 end of out	

    EE.VLD.128.IP q0, a2, 16
    EE.VLD.128.IP q1, a2, 16

    .big_loop_start:
        EE.VUNZIP.16 q0, q1         // x -> q0, y -> q1
        EE.ZERO.Q q6            // sum -> q6
        EE.MOV.S16.QACC q0      // x-> QACC
        EE.VRELU.S16 q0, a6, a7     // abs(x) -> q0
        MOVI.N a8, 0
        loopnez a9, .loop_cordic_end
            // step 0
            EE.VCMP.LT.S16 q5, q1, q7  // neg y
            EE.VCMP.GT.S16 q4, q1, q7  // pos y
            // xn = x1 + (y >> a8) or xn = x1 - (y >> a8)
            WSR.SAR a8
            EE.VMUL.S16 q3, q1, q4 // pos -y
            EE.VSUBS.S16 q2, q0, q3 // |x|-(-y) -> q2
            EE.VMUL.S16 q3, q1, q5 // neg -y
            EE.VADDS.S16 q2, q2, q3 // |x|+(-y) -> q2
            // y = y - (x1 >> a8) or y = y + (x1 >> a8)
            EE.VMUL.S16 q3, q0, q4 // pos -|x|
            EE.VADDS.S16 q1, q1, q3 // y = y + (-|x|); -> q1
            EE.VMUL.S16 q3, q0, q5 // neg -|x|
            EE.VSUBS.S16 q1, q1, q3 // y - (-|x|) -> q1
            // x1 = xn
            EE.VADDS.S16 q0, q2, q7     
            //  sum += tan[a8] or sum -= tan[a8]
            WSR.SAR a7
            EE.VLDBC.16.IP q2, a5, 2    // q2 = 8{tan[a8]}
            EE.VMUL.S16 q4, q2, q4      // pos -tan[a8]
            EE.VSUBS.S16 q6, q6, q4     // sum = sum -(-tan[a8])
            EE.VMUL.S16 q5, q2, q5      // neg -tan[a8]
            EE.VADDS.S16 q6, q6, q5     // sum = sum +(-tan[a8])

            ADDI.N a8,a8,1
        .loop_cordic_end:
        
        EE.SRCMB.S16.QACC q3, a7, 0 // x        -> q3
        EE.VLDBC.16.IP q2, a5, 2    // -16383   -> q2
        EE.VLDBC.16.IP q1, a5, 2    // 1        -> q1
        EE.VCMP.LT.S16 q5, q3, q7   // neg x     -> q5
        EE.VCMP.LT.S16 q0, q6, q7   // neg sum   -> q0
        EE.VMUL.S16 q4, q5, q0      // (neg x) & (neg sum) -> q4
        EE.VMUL.S16 q4, q2, q4      // erase and (-16383)

        EE.VADDS.S16 q0, q0, q1     // not neg sum   -> q0
        EE.VMUL.S16 q0, q5, q0      // (neg x) & (not neg sum) -> q4
        EE.VMUL.S16 q0, q2, q0      // erase and (16383)
        EE.VADDS.S16 q3, q0, q4     // combine (-16383) and (16383) -> q3

        EE.VCMP.EQ.S16 q4, q3, q7   // not neg x    -> q4
        EE.VLDBC.16 q5, a5          // 2            -> q5
        EE.VMUL.S16.LD.INCP q0, a2, q4, q4, q5
        EE.VADDS.S16 q4, q4, q1     // sign for sum
        EE.VMUL.S16.LD.INCP q1, a2, q6, q6, q4

        EE.VSUBS.S16 q6, q3, q6     // (-16383) and (16383) - sum

        EE.VST.128.IP q6, a3, 16 //save arqument
        MOV.N a5, a10

        BLT a3, a4, .big_loop_start

    retw.n
