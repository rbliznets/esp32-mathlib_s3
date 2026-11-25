    .section .iram1,"ax"  // Place this function in the IRAM1 section (for faster access)
    .text
    .align  4
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Rotation of a complex q15 vector by a periodic twiddle table:
//     data[i] = data[i] * tw[(phase + i) % period]
// EE.CMUL.S16 processes one 64-bit half of the operands per instruction (immediate 0 - complex 0, 1;
// 1 - complex 2, 3, the other half of the destination is kept), so two instructions rotate 4 complex
// samples. The table is stored per phase - tw[p] holds the four twiddles of the phases p, p+1, p+2,
// p+3 (modulo period), so a 128-bit load at tw[phase] always gives four consecutive twiddles even
// when the period is not a multiple of 4 (the modem rotates by -pi/3 per sample at 48 kHz, period 6).
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   rotate_q15_pie,@function
    .align  4
    .global rotate_q15_pie
rotate_q15_pie:
// data   - a2  (complex_q15 vector, 16 bytes aligned, modified in place)
// size   - a3  (vector size, multiple of 4)
// tw     - a4  (phase table, 16 bytes aligned: tw[p][j] = twiddle of phase (p + j) % period)
// period - a5  (number of phases, >= 1)
// phase  - a6  (starting phase, 0..period-1)
    entry   a1, 16

    ssai    15                      // Q15 scaling for EE.CMUL.S16
    srli    a3, a3, 2               // number of groups of 4 complex samples
    slli    a6, a6, 4               // phase -> byte offset in the table
    slli    a5, a5, 4               // period -> table size in bytes
    mov.n   a9, a2                  // store pointer
    movi.n  a10, 64                 // 4 phases per iteration

    beqz    a3, .Lend_rotate_q15_pie
    beqz    a5, .Lend_rotate_q15_pie
    bltu    a10, a5, .Lloop_rotate_q15_pie
    remu    a10, a10, a5            // (4 % period) rows; one subtraction now suffices below
.Lloop_rotate_q15_pie:
    add.n   a7, a4, a6
    ee.vld.128.ip   q1, a7, 0       // twiddles of the phases phase..phase+3
    ee.vld.128.ip   q0, a2, 16      // 4 complex samples
    ee.cmul.s16     q2, q0, q1, 0   // complex 0, 1 (lanes 0..1 of the destination)
    ee.cmul.s16     q2, q0, q1, 1   // complex 2, 3 (keeps lanes 0..1 of q2)
    ee.vst.128.ip   q2, a9, 16
    add.n   a6, a6, a10
    bltu    a6, a5, .Lwrap_done_rotate_q15_pie
    sub     a6, a6, a5              // wrap the phase offset around the table
.Lwrap_done_rotate_q15_pie:
    addi.n  a3, a3, -1
    bnez    a3, .Lloop_rotate_q15_pie

.Lend_rotate_q15_pie:
    retw.n  // Return from function (windowed return)
