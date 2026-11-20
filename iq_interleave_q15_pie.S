    .section .iram1,"ax"  // Place these functions in the IRAM1 section (for faster access)
    .text
    .align  4
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Split and merge of interleaved complex q15 vectors (I, Q, I, Q, ...).
// EE.VUNZIP.16 puts the even 16-bit lanes of both source vectors into the first destination and the
// odd lanes into the second one, EE.VZIP.16 does the reverse, so both directions handle 8 complex
// samples per iteration.
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   iq_unzip_q15_pie,@function
    .align  4
    .global iq_unzip_q15_pie
iq_unzip_q15_pie:
// in   - a2  (complex_q15 vector, 16 bytes aligned)
// re   - a3  (output q15 vector, 16 bytes aligned)
// im   - a4  (output q15 vector, 16 bytes aligned)
// size - a5  (vector size, multiple of 8)
    entry   a1, 16

    srli    a5, a5, 3               // number of groups of 8 complex samples
    loopnez a5, .Lend_iq_unzip_q15_pie
        ee.vld.128.ip   q0, a2, 16  // re0, im0, ... re3, im3
        ee.vld.128.ip   q1, a2, 16  // re4, im4, ... re7, im7
        ee.vunzip.16    q0, q1      // q0 = re[0..7], q1 = im[0..7]
        ee.vst.128.ip   q0, a3, 16
        ee.vst.128.ip   q1, a4, 16
    .Lend_iq_unzip_q15_pie:

    retw.n  // Return from function (windowed return)
///////////////////////////////////////////////////////////////////////////////////////////////////////
    .type   iq_zip_q15_pie,@function
    .align  4
    .global iq_zip_q15_pie
iq_zip_q15_pie:
// re   - a2  (input q15 vector, 16 bytes aligned)
// im   - a3  (input q15 vector, 16 bytes aligned)
// out  - a4  (output complex_q15 vector, 16 bytes aligned)
// size - a5  (vector size, multiple of 8)
    entry   a1, 16

    srli    a5, a5, 3               // number of groups of 8 complex samples
    loopnez a5, .Lend_iq_zip_q15_pie
        ee.vld.128.ip   q0, a2, 16  // re[0..7]
        ee.vld.128.ip   q1, a3, 16  // im[0..7]
        ee.vzip.16      q0, q1      // q0 = re0, im0, ... re3, im3; q1 = re4, im4, ... re7, im7
        ee.vst.128.ip   q0, a4, 16
        ee.vst.128.ip   q1, a4, 16
    .Lend_iq_zip_q15_pie:

    retw.n  // Return from function (windowed return)
