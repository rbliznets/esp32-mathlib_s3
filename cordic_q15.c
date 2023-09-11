#include <assert.h>
#include "mathlib_s3.h"
#include "complex_s3.h"
#include "sdkconfig.h"

//static const q15 tan_array[16]={4095, 2417, 1277, 648, 325, 162, 81, 40, 20, 10, 4, 2, 1, -16383, 1, 2};
static const q15 tan_array[16]={4096, 2418, 1277, 648, 325, 163, 81, 41, 20, 10, 5, 2, 1, -16383, 1, 2};
q15 atan2_q15_s3(q15 y, q15 x, const q15* tan);
inline q15 atan2_q15(q15 y, q15 x)
{
    return atan2_q15_s3(y,x,tan_array);
}

void arg_16_q15_pie(complex_q15* in, q15* out, uint32_t size, const q15* tan);
inline void arg_16_q15(complex_q15* in, q15* out, uint32_t size)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)in % 16) == 0);
    assert(((uint32_t)out % 16) == 0);
    assert((size % 8) == 0);
    assert(size > 0);
#endif 
    arg_16_q15_pie(in, out, size, tan_array);
}

void sincos_q15_s3(q15 angle, q15* cs, q15* sn, const q15* tan);
inline void sincos_q15(q15 angle, q15* sn, q15* cs)
{
    sincos_q15_s3(angle, cs, sn, tan_array);
}
