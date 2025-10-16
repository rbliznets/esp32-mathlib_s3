    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access)
    .text
    .type   cmul_q15_pie,@function
    .align  4
    .global cmul_q15_pie
cmul_q15_pie: 
// x   - a2  (first complex number input)
// y   - a3  (second complex number input)
// return z=x*y - a2 (result returned in a2)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI
    
    ssai    15      // Set Shift Amount for Immediate instructions to 15

    EE.MOVI.32.Q q0, a2, 0  // Move the 32-bit complex value from a2 into vector register q0
    EE.MOVI.32.Q q1, a3, 0  // Move the 32-bit complex value from a3 into vector register q1
    EE.CMUL.S16 q2, q0, q1, 0  // Perform complex multiplication: q2 = q0 * q1 (S16 signed 16-bit)
    EE.MOVI.32.A q2, a2, 0     // Move the 32-bit complex result from q2 back to register a2 for return

    retw.n  // Return from function (windowed return)
/////////////////////////////////////////////////////////////////////////////////
    .text
    .type   cmul10_q15_pie,@function
    .align  4
    .global cmul10_q15_pie
cmul10_q15_pie: 
// in1[10]   - a2  (pointer to input vector of 10 complex numbers)
// k[1]      - a3  (pointer to scalar complex number)
// out[10]   - a4  (pointer to output vector of 10 complex numbers)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI
    
    ssai    15              // Set Shift Amount for Immediate instructions to 15
    EE.VLDBC.32 q1, a3      // Broadcast the 32-bit scalar complex number from a3 into vector register q1 (so q1 contains 4 copies of k)

    EE.VLD.128.IP q0, a2, 16    // Load 128 bits (4 complex numbers) from in1[0..3] into q0, increment a2 by 16
    EE.CMUL.S16 q2, q0, q1, 0   // Multiply first 4 complex numbers: q2 = q0 * q1
    EE.CMUL.S16.LD.INCP q3, a2, q2, q0, q1, 1 // Multiply q0 * q1 -> q2, load next 4 complex numbers in1[4..7] -> q3, increment a2 by 16

    EE.CMUL.S16.ST.INCP q2, a4, q4, q3, q1, 0 // Multiply q3 * q1 -> q4, store previous result q2 to out[0..3] -> a4, increment a4 by 16
    EE.CMUL.S16.LD.INCP q0, a2, q4, q3, q1, 1 // Multiply q3 * q1 -> q4, load next 4 complex numbers in1[8..11] -> q0, increment a2 by 16
 
    EE.CMUL.S16.ST.INCP q4, a4, q2, q0, q1, 0 // Multiply q0 * q1 -> q2, store previous result q4 to out[4..7] -> a4, increment a4 by 16
    EE.VST.L.64.IP q2, a4, 0   // Store the low 64 bits (last 2 complex numbers) from q2 to out[8..9], increment a4 by 8

    retw.n  // Return from function (windowed return)
