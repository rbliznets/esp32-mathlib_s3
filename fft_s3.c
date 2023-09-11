/*!
	\file
	\brief Оптимизированные функции для FFT.
	\authors Близнец Р.А.
	\version 0.0.0.1
	\date 13.03.2023
*/

#include "fft_s3.h"
#include <math.h>

#include <stdio.h>

int16_t fft_log2(uint32_t fftSize)
{
    uint32_t x = 0x8000;
    int16_t res = 15;
    while (x != 0) 
    {
        if( ((x & fftSize) != 0) || (res==0))
        {
            if (x == fftSize)
            {
                return res;
            } 
            else
            {
                return -res;
            }
        }
        res--;
        x >>= 1;
    }
    return 0;
}

inline uint32_t revbin_update(uint32_t r, uint32_t n)
{
    for (uint32_t m=n>>1; (!((r^=m)&m)); m>>=1);
    return r;
}

void revbin_permute(complex_q15* data, uint32_t fftSize)
{
    uint32_t* dt = (uint32_t*)data;
    uint32_t nh = fftSize >> 1;
    uint32_t r = 0;
    uint32_t x = 1;
    uint32_t t;
    while( x < nh)
    {
        r = r + nh;;
        t = dt[x];
        dt[x] = dt[r];
        dt[r] = t;
        x++;

        r = revbin_update(r,fftSize);
        if (r > x)
        {
            t = dt[x];
            dt[x] = dt[r];
            dt[r] = t;
            t = dt[fftSize-1-x];
            dt[fftSize-1-x] = dt[fftSize-1-r];
            dt[fftSize-1-r] = t;
        }
        x++;
   }
}

// fftSize - 2
void init_fft(complex_q15* w, uint32_t fftSize)
{
    int16_t n = fft_log2(fftSize);
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)w % 16) == 0);
    assert(n >= 3);
#endif 
    float e = M_PI * 2.0 / fftSize;
    fftSize >>= 1;
    for (int i = 0; i < fftSize; i++) 
    {
        w[i].re = (q15)roundf(INT16_MAX * cosf(i * e));
        w[i].im = (q15)roundf(INT16_MAX * sinf(i * e));
    }

    complex_q15* w_last = w;
    complex_q15* w_cur = &w_last[fftSize];
    for (int i = 3; i < n; i++) 
    {
        fftSize >>= 1;
        for(int j = 0; j < fftSize; j++)
        {
           w_cur[j].re = w_last[2*j].re; 
           w_cur[j].im = w_last[2*j].im; 
        }
        w_last = w_cur;
        w_cur = &w_last[fftSize];
    }
    fftSize >>= 1;
    for(int j = 0; j < fftSize; j++)
    {
        w_cur[j].re = -w_last[2*j].re; 
        w_cur[j].im = -w_last[2*j].im; 
    }
}

complex_q15* getW(complex_q15* w, uint32_t fftSize, uint32_t fftSize2)
{
#ifdef CONFIG_CHECK_PARAM
    assert(fftSize >= fftSize2) == 0);
#endif
    complex_q15* w_cur = w;
    while(fftSize > fftSize2)
    {
        fftSize >>= 1;
        w_cur = &w_cur[fftSize];
    } 
    return w_cur;
}

void fft_r2_q15_pie(complex_q15* data, complex_q15* w, uint32_t fftSize);
inline void fft_radix2(complex_q15* data, complex_q15* w, uint32_t fftSize)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)data % 16) == 0);
    assert(fftSize >= 16);
#endif 
    fft_r2_q15_pie(data, w, fftSize);
}
 

void fft_r2s_q15_pie(complex_q15* data, complex_q15* w, uint32_t fftSize);
inline void fft_radix2_scale(complex_q15* data, complex_q15* w, uint32_t fftSize)
{
#ifdef CONFIG_CHECK_PARAM
    assert(((uint32_t)data % 16) == 0);
    assert(fftSize >= 16);
#endif 
    fft_r2s_q15_pie(data, w, fftSize);
}
