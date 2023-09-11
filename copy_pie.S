    .text
    .type   copy_pie,@function
    .align  4
    .global copy_pie
copy_pie: 
// in   - a2
// out   - a3
// size  - a4
    entry	a1, 16

    srli 	a4, a4, 3           //  size/8 	

    EE.LD.128.USAR.IP   q0, a2, 16 
    
    loopnez a4, .loop_end_copy_pie
        EE.LD.128.USAR.IP   q1, a2, 16 
        EE.SRC.Q.QUP        q2, q0,q1
        ee.vst.128.ip       q2, a3, 16
    .loop_end_copy_pie:

    retw.n
/////////////////////////////////////////////////////////////////////////////////
    .text
    .type   copy_16_pie,@function
    .align  4
    .global copy_16_pie
copy_16_pie: 
// in   - a2
// out   - a3
// size  - a4
    entry	a1, 16

    srli 	a4, a4, 3           //  size/8 	

    loopnez a4, .loop_end_copy_16_pie
        ee.vld.128.ip   q1, a2, 16 
        ee.vst.128.ip   q1, a3, 16
    .loop_end_copy_16_pie:

    retw.n
