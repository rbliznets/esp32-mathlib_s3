/*!
	\file
	\brief Оптимизированные функции DSP.
	\authors Близнец Р.А.
	\version 0.3.0.0
	\date 11.11.2022
*/

#include <assert.h>
#include "mathlib_s3.h"
#include "sdkconfig.h"

void copy_pie(q15* in, q15* out, uint32_t size);
inline void copy(q15* in, q15* out, uint32_t size)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)out % 16) == 0);
    assert((size % 8) == 0);
    assert(size > 0);
#endif 
    copy_pie(in, out, size);
}

void copy_16_pie(q15* in, q15* out, uint32_t size);
inline void copy_16(q15* in, q15* out, uint32_t size)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)in % 16) == 0);
    assert(((uint32_t)out % 16) == 0);
    assert((size % 8) == 0);
    assert(size > 0);
#endif 
    copy_16_pie(in, out, size);
}

void scaleVector_q15_pie(q15* in, q15* k, q15* out, uint32_t size);
inline void scaleVector(q15* in, q15* k, q15* out, uint32_t size)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)k % 2) == 0);
    assert(((uint32_t)out % 16) == 0);
    assert((size % 8) == 0);
    assert(size >= 16);
#endif 
    scaleVector_q15_pie(in, k, out, size);
}

void scaleVector_q15_16_pie(q15* in, q15* k, q15* out, uint32_t size);
inline void scaleVector_16(q15* in, q15* k, q15* out, uint32_t size)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)in % 16) == 0);
    assert(((uint32_t)k % 2) == 0);
    assert(((uint32_t)out % 16) == 0);
    assert((size % 8) == 0);
    assert(size > 0);
#endif 
    scaleVector_q15_16_pie(in, k, out, size);
}

void shrinkVector_16_pie(uint32_t* in, uint8_t shift, q15* out, uint32_t size);
inline void shrinkVector_16(uint32_t* in, uint8_t shift, q15* out, uint32_t size)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)in % 16) == 0);
    assert(shift < 32 );
    assert(((uint32_t)out % 16) == 0);
    assert((size % 8) == 0);
    assert(size > 0);
#endif 
    shrinkVector_16_pie(in, shift, out, size);
}

q15 dot_product_q15_16_16_pie(q15* in1, q15* in2, uint32_t size);
inline q15 dot_product_16_16(q15* in1, q15* in2, uint32_t size)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)in1 % 16) == 0);
    assert(((uint32_t)in2 % 16) == 0);
    assert((size % 8) == 0);
    assert(size > 0);
#endif 
    return dot_product_q15_16_16_pie(in1, in2, size);
}

q15 dot_product_q15_1_16_pie(q15* in1, q15* in2, uint32_t size);
inline q15 dot_product_1_16(q15* in1, q15* in2, uint32_t size)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)in2 % 16) == 0);
    assert((size % 8) == 0);
    assert(size > 0);
#endif 
    return dot_product_q15_1_16_pie(in1, in2, size);
}

void fir_1_16(q15* in, q15* k, uint32_t ksize, q15* out, uint32_t size)
{
#ifdef CONFIG_CHECK_PARAM
    assert(size > 0);
#endif 
    for(uint32_t i=0; i < size; i++) out[i] = dot_product_1_16(&in[i],k,ksize);
}

void addVectors_q15_pie(q15* in1, q15* in2, q15* out, uint32_t size);
inline void addVectors_q15(q15* in1, q15* in2, q15* out, uint32_t size)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)in1 % 16) == 0);
    assert(((uint32_t)in2 % 16) == 0);
    assert(((uint32_t)out % 16) == 0);
    assert((size % 8) == 0);
    assert(size > 0);
#endif 
    addVectors_q15_pie(in1, in2, out, size);
}

int16_t normalize_q15_pie(q15* in, q15* out, uint32_t size);
inline int16_t normalize_q15(q15* in, q15* out, uint32_t size)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)in % 16) == 0);
    assert(((uint32_t)out % 16) == 0);
    assert((size % 8) == 0);
    assert(size > 0);
#endif 
    return normalize_q15_pie(in, out, size);
}

int16_t normalize_q14_pie(q15* in, q15* out, uint32_t size);
inline int16_t normalize_q14(q15* in, q15* out, uint32_t size)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)in % 16) == 0);
    assert(((uint32_t)out % 16) == 0);
    assert((size % 8) == 0);
    assert(size > 0);
#endif 
    return normalize_q14_pie(in, out, size);
}

void fir_16_16_q15_pie(q15* in, q15* k, uint32_t ksize, q15* out, uint32_t size);
inline void fir_16_16_q15(q15* in, q15* k, uint32_t ksize, q15* out, uint32_t size)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)in % 16) == 0);
    assert(((uint32_t)k % 16) == 0);
    assert(((uint32_t)out % 16) == 0);
    assert((ksize % 8) == 0);
    assert(ksize > 0);
    assert((size % 8) == 0);
    assert(size > 8);
#endif 
    fir_16_16_q15_pie(in, k, ksize, out, size);
}
