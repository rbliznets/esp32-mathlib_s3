/*!
	\file
	\brief Оптимизированные функции для комплексных чисел.
	\authors Близнец Р.А.
	\version 0.3.0.0
	\date 14.02.2023
	\copyright (c) Copyright 2023, ООО "Глобал Ориент", Москва, Россия, https://github.com/Global-Orient/Gonets
*/

#include <assert.h>
#include "complex_s3.h"
#include "sdkconfig.h"


void magnitude_q15_pie(complex_q15* in, q15* out, uint32_t size);
inline void magnitude_q15(complex_q15* in, q15* out, uint32_t size)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)in % 16) == 0);
    assert(((uint32_t)out % 16) == 0);
    assert((size % 8) == 0);
    assert(size > 0);
#endif 
    magnitude_q15_pie(in, out, size);
}

complex_q15 cmul_q15_pie(complex_q15 x, complex_q15 y);
inline complex_q15 cmul_q15(complex_q15 x, complex_q15 y)
{
	return cmul_q15_pie(x, y);
};

void cmul10_q15_pie(complex_q15* in, complex_q15* k, complex_q15* out);
void cmul10_q15(complex_q15* in, complex_q15* k, complex_q15* out)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)in % 16) == 0);
    assert(((uint32_t)k % 4) == 0);
    assert(((uint32_t)out % 16) == 0);
#endif 
    cmul10_q15_pie(in, k, out);
}
