    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access)
    .text
    .type   copy_pie,@function
    .align  4
    .global copy_pie
copy_pie: 
// in   - a2  (pointer to input array)
// out  - a3  (pointer to output array)
// size - a4  (number of elements to copy)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    srli 	a4, a4, 3           // Divide size by 8 (as the loop processes 8 elements per iteration)

    EE.LD.128.USAR.IP   q0, a2, 16  // Load 128 bits (8 int16_t) from input (a2) using SAR shift, increment a2 by 16
    
    loopnez a4, .loop_end_copy_pie  // Loop 'size/8' times, exit if a4 becomes zero
        EE.LD.128.USAR.IP   q1, a2, 16  // Load next 128 bits (8 int16_t) from input (a2) using SAR shift, increment a2 by 16
        EE.SRC.Q.QUP        q2, q0,q1   // Shift and combine q0 and q1 (likely interleaving or reordering data into q2)
        ee.vst.128.ip       q2, a3, 16  // Store 128 bits (8 int16_t) from q2 to output (a3), increment a3 by 16
    .loop_end_copy_pie:

    retw.n  // Return from function (windowed return)
/////////////////////////////////////////////////////////////////////////////////
    .text
    .type   copy_16_pie,@function
    .align  4
    .global copy_16_pie
copy_16_pie: 
// in   - a2  (pointer to input array, 16-byte aligned)
// out  - a3  (pointer to output array, 16-byte aligned)
// size - a4  (number of elements to copy)
    entry	a1, 16  // Function entry, allocate 16 bytes on the stack for ABI

    srli 	a4, a4, 3           // Divide size by 8 (as the loop processes 8 elements per iteration)

    loopnez a4, .loop_end_copy_16_pie  // Loop 'size/8' times, exit if a4 becomes zero
        ee.vld.128.ip   q1, a2, 16  // Load 128 bits (8 int16_t) from input (a2), increment a2 by 16
        ee.vst.128.ip   q1, a3, 16  // Store 128 bits (8 int16_t) from q1 to output (a3), increment a3 by 16
    .loop_end_copy_16_pie:

    retw.n  // Return from function (windowed return)
