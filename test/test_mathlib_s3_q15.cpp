/**
 * @file test_mathlib_s3_q15.cpp
 * @author Roman Bliznets (r.bliznets@gmail.com)
 * @brief Unity module tests with benchmark
 * @version 0.0.0.1
 * @date 22.09.2022
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <limits.h>
#include <cstring>
#include "unity.h"
#include "mathlib_s3.h"
#include "complex_s3.h"
#include "fft_s3.h"
#include "CTrace.h"
#include <math.h>
#include "esp_dsp.h"

#define countof(x) (sizeof(x)/sizeof(x[0]))
#define N   10000 /**<  Number of iterations */
#define TEST_SIZE (1024)

#define FFT_SIZE 4096
__attribute__((aligned(16)))
static complex_q15 w[FFT_SIZE];//FFT_SIZE-4

__attribute__((aligned(16)))
static q15 in1[TEST_SIZE];
__attribute__((aligned(16)))
static q15 in2[TEST_SIZE];
__attribute__((aligned(16)))
static complex_q15 in[FFT_SIZE+16];
__attribute__((aligned(16)))
static q15 out[TEST_SIZE];

TEST_CASE("fft_q15", "[math][mathlib_s3]")
{
   TEST_ASSERT_EQUAL_INT16(10,fft_log2(1024));
   TEST_ASSERT_EQUAL_INT16(-9,fft_log2(1000));
   TEST_ASSERT_EQUAL_INT16(0,fft_log2(0));

   std::memset(w,0,sizeof(w));
   STARTTIMESHOT();
   init_fft(w,FFT_SIZE);
   STOPTIMESHOT("init_fft time");
   //TRACEDATA("w",(int16_t*)(&w[FFT_SIZE-8]),16);

   STARTTIMESHOT();
   dsps_bit_rev_sc16_ansi((int16_t*)in,FFT_SIZE);
   STOPTIMESHOT("dsps_bit_rev_sc16_ansi time");
   STARTTIMESHOT();
   revbin_permute(in,FFT_SIZE);
   STOPTIMESHOT("revbin_permute time");

   float ex = M_PI * 2.0 / 10;
   std::memset(in,0,sizeof(in));
   for(int i = 0; i < FFT_SIZE; i++)
   {
      in[i].re = (q15)(INT16_MAX/2 * sinf(i * ex));
      in[i].im = (q15)(INT16_MAX/2 * cosf(i * ex));
   }
   TDEC("FFT_SIZE",FFT_SIZE);
   dsps_fft2r_init_sc16(nullptr, FFT_SIZE);
   STARTTIMESHOT();
   //dsps_fft2r_sc16((int16_t*)in, FFT_SIZE);
   fft_radix2(in, w, FFT_SIZE);
   // fft_radix2_scale(in, w, FFT_SIZE);
   STOPTIMESHOT("fft_radix2_scale time");
   //TRACEDATA("in0.7",(int16_t*)in,8);
   //TRACEDATA("in8.15",(int16_t*)&in[4],8);
   //TRACEDATA("in16",(int32_t*)&in[0],1);

   revbin_permute(in,FFT_SIZE);
   magnitude_q15(&in[FFT_SIZE-(TEST_SIZE/2)],&out[0],TEST_SIZE/2);
   magnitude_q15(&in[0],&out[TEST_SIZE/2],TEST_SIZE/2);
   TRACEDATA("out",out,TEST_SIZE);
}


/**
 * @brief Sum of 2 vectors
 * 
 * @param in1 1 vector
 * @param in2 2 vector
 * @param out sum
 * @param size size of vector
 */
static void addVectors(q15* in1, q15* in2, q15* out, uint32_t size)
{
   for(uint32_t i=0; i<size; i++)
   {
      out[i]=in1[i]+in2[i];
   }
}

TEST_CASE("addVectors_q15", "[math][mathlib_s3]")
{
   uint32_t mem1=esp_get_free_heap_size();

   TRACE("Number of iterations",N,false);
   TRACE("Vector size",countof(in1),false);
   for(int16_t i=0; i < countof(in1); i++)
   {
      in1[i]=i;
      in2[i]=i;
   }
   
   STARTTIMESHOT();
   for(int i=0;i<N;i++)
   {
      addVectors(in1,in2,out,countof(in1));
   }
   STOPTIME("C++ version time",N);

   for( auto&x : out) x=0;
   STARTTIMESHOT();
   for(int i=0;i<N;i++)
   {
      addVectors_q15(in1,in2,out,countof(in1));
   }
   STOPTIME("addVectors_q15 time",N);

   for(int16_t i=0; i < countof(in1); i++)
   {
      TEST_ASSERT_EQUAL_INT16(i+i, out[i]);
   }

   for( auto&x : out) x=0;
   addVectors_q15(in1,in2,out,8);
   for(int16_t i=0; i < 8; i++)
   {
      TEST_ASSERT_EQUAL_INT16(i+i, out[i]);
   }

   uint32_t mem2=esp_get_free_heap_size();
   if(mem1 != mem2)
   {
      TRACE("memory leak",mem1-mem2,false);
      TRACE("start",mem1,false);
      TRACE("stop",mem2,false);
      TEST_FAIL_MESSAGE("memory leak");
}
}

static void magnitude(complex_q15* in, q15* out, uint32_t size)
{
    for(uint32_t i=0; i<size; i++)
    {
        int32_t x=in[i].re;
        int32_t y=in[i].im;
        out[i]=(((x*x) >> 16) + ((y*y) >> 16));
    }
}

TEST_CASE("magnitude_q15", "[math][mathlib_s3]")
{
   uint32_t mem1=esp_get_free_heap_size();

   TRACE("Number of iterations",N,false);
   TRACE("Input size",countof(in),false);
   for(int16_t i=0; i < countof(in); i++)
   {
        in[i].re = INT16_MAX/(i+1);
        in[i].im = INT16_MAX/(i+1);
   }
   in[countof(in)-1].re = INT16_MAX;
   in[countof(in)-1].im = INT16_MAX;
   
   STARTTIMESHOT();
   magnitude(in,in1,countof(in));
   STOPTIME("C++ version 1 time",1);
   STARTTIMESHOT();
   for(int i=0;i<N;i++)
   {
      magnitude(in,in1,countof(in));
   }
   STOPTIME("C++ version time",N);

   STARTTIMESHOT();
   magnitude_q15(in,out,countof(in));
   STOPTIME("magnitude_q15 1 time",1);
   for( auto& x : out) x=0;
   STARTTIMESHOT();
   for(int i=0;i<N;i++)
   {
      magnitude_q15(in,out,countof(in));
   }
   STOPTIME("magnitude_q15 time",N);

   TEST_ASSERT_EQUAL_INT16_ARRAY(out, in1, countof(in));

   for( auto& x : out) x=0;
   magnitude_q15(in,out,8);
   TEST_ASSERT_EQUAL_INT16_ARRAY(out, in1, 8);

   uint32_t mem2=esp_get_free_heap_size();
   if(mem1 != mem2)
   {
      TRACE("memory leak",mem1-mem2,false);
      TRACE("start",mem1,false);
      TRACE("stop",mem2,false);
      TEST_FAIL_MESSAGE("memory leak");
}
}

static int16_t normalize(q15* _in, q15* _out, uint32_t size)
{
    q15 mn=0;
    q15 mx=0;
    for(uint32_t i=0; i<size; i++ )
    {
        if(_in[i] < mn) mn = _in[i];
        if(_in[i] > mx) mx = _in[i];
    }

    if(mx < (-mn)) mx=-mn;
    if(mx == 0) return 0;
    int16_t mlt= (0x7fff/mx);
    if(mlt == 1) return 1;

    for(uint32_t i=0; i<size; i++ )
    {
        _out[i]=_in[i]*mlt;
    }
    return mlt;
}

TEST_CASE("normalize_q15", "[math][mathlib_s3]")
{
   uint32_t mem1=esp_get_free_heap_size();

   TRACE("Number of iterations",N,false);
   TRACE("Input size",countof(in1),false);
   for(int16_t i=0; i < countof(in1); i++)
   {
        in1[i] = INT16_MAX/(i+3);
   }
   in[countof(in1)-1].re = INT16_MAX/3;
   
   int16_t mult=0;
   STARTTIMESHOT();
   mult=normalize(in1,in2,countof(in1));
   STOPTIME("C++ version 1 time",1);
   // TRACE("mult",mult,false);

   STARTTIMESHOT();
   for(int i=0;i<N;i++)
   {
      mult=normalize(in1,in2,countof(in1));
   }
   STOPTIME("C++ version time",N);
   // TRACEDATA("in2",in2,8);
   
   std::memset(out,0,sizeof(out));
   TRACEDATA("in1",in1,8);
   STARTTIMESHOT();
   mult=normalize_q15(in1,out,countof(in1));
   STOPTIME("normalize_q15 1 time",1);
   TRACEDATA("out",out,8);

   for( auto&x : out) x=0;
   STARTTIMESHOT();
   for(int i=0;i<N;i++)
   {
      mult=normalize_q15(in1,out,countof(in1));
   }
   STOPTIME("normalize_q15 time",N);
   TRACE("mult",mult,false);
   // TRACEDATA("out",out,8);

   TEST_ASSERT_EQUAL_INT16_ARRAY(out, in2, countof(in2));

   uint32_t mem2=esp_get_free_heap_size();
   if(mem1 != mem2)
   {
      TRACE("memory leak",mem1-mem2,false);
      TRACE("start",mem1,false);
      TRACE("stop",mem2,false);
      TEST_FAIL_MESSAGE("memory leak");
   }
}

inline q15 dot_product_c(q15* din1, q15* din2, uint32_t size)
{
    int32_t sum=0;
    for (uint32_t i = 0; i < size; i++)
    {
        sum += ((int32_t)din1[i])*din2[i];
    }
    return (q15)(sum >> 15);
}

TEST_CASE("dot_product", "[math][mathlib_s3]")
{
   uint32_t mem1=esp_get_free_heap_size();

   TDEC("**** Data size ****",TEST_SIZE-8);
   TDEC("Number of iterations",N);

   for(int16_t i=0; i < TEST_SIZE; i++)
   {
      in1[i]=toQ15(0.5)/16;
      in2[i]=toQ15(0.03);
      out[i]=0;
   }
   in1[1]=toQ15(-0.25);

   STARTTIMESHOT();
   out[0]=dot_product_16_16(in1, in2, TEST_SIZE-8);
   STOPTIME("dot_product_16_16 1 time",1);
   TEST_ASSERT_EQUAL_INT16(out[0], 30903);
   STARTTIMESHOT();
   for(int i=0; i<N; i++) out[0]=dot_product_16_16(in1, in2, TEST_SIZE-8);
   STOPTIME("dot_product_16_16 time",N);

   STARTTIMESHOT();
   out[1]=dot_product_1_16(&in2[1], in2, TEST_SIZE-8);
   STOPTIME("dot_product_1_16 1 time",1);
   TEST_ASSERT_EQUAL_INT16(out[1], 29960);
   STARTTIMESHOT();
   for(int i=0; i<N; i++) out[1]=dot_product_1_16(&in2[1], in2, TEST_SIZE-8);
   STOPTIME("dot_product_1_16 time",N);

   STARTTIMESHOT();
   out[2]=dot_product_c(in1, in2, TEST_SIZE-8);
   STOPTIME("dot_product_c 1 time",1);
   STARTTIMESHOT();
   for(int i=0; i<N; i++) out[2]=dot_product_c(in1, in2, TEST_SIZE-8);
   STOPTIME("dot_product_c time",N);

   uint32_t mem2=esp_get_free_heap_size();
   if(mem1 != mem2)
   {
      TRACE("memory leak",mem1-mem2,false);
      TRACE("start",mem1,false);
      TRACE("stop",mem2,false);
      TEST_FAIL_MESSAGE("memory leak");
   }
}

#include "dsps_dotprod.h"
__attribute__((aligned(16)))
static const int16_t gauss[24]={1,4,15,52,151,380,832,1579,2599,3710,4594,4933,4594,3710,2599,1579,832,380,151,52,15,4,1,0};
void fir_dsp(int16_t* inData, int16_t* outData, uint16_t data_size, const int16_t* k, uint16_t k_size)
{
   for(int i=0; i < data_size; i++)
   {
      dsps_dotprod_s16_ae32(&inData[i], k, &outData[i], k_size, 0);
   }
}
void fir_q15_asc(int16_t* inData, int16_t* outData, uint16_t data_size, const int16_t* k, uint16_t k_size)
{
    int32_t sum;
    for(int i=0; i < data_size; i++)
    {
        sum=0;
        for(int j=0; j < k_size; j++)
        {
            sum += inData[i+j] * k[j];
        }
        outData[i] = (int16_t)(sum >> 15);
    }
}

TEST_CASE("fir_q15", "[math][mathlib_s3]")
{
   uint32_t mem1=esp_get_free_heap_size();

   // TDEC("Number of iterations",N);
   uint16_t size = (sizeof(in1)/(sizeof(int16_t)))-24;
   TDEC("***** Data size *******",size);
   int16_t* dataIn = (int16_t*)in1;
   int16_t* dataOut = (int16_t*)out;
   for(int i=0;i<(size/20);i++)
   {
      for (size_t j = 0; j < 10; j++)
      {
         dataIn[i*20+j] = 10000;
      }
      for (size_t j = 10; j < 20; j++)
      {
         dataIn[i*20+j] = -10000;
      }
   }
   std::memset(out,0,sizeof(out));

   STARTTIMESHOT();
   fir_q15_asc(dataIn, dataOut, size, gauss, 23);
   STOPTIME("fir_q15_asc", 1);
   STARTTIMESHOT();
   for(int i=0; i<N; i++) fir_q15_asc(dataIn, dataOut, size, gauss, 23);
   STOPTIME("fir_q15_asc time",N);

   std::memset(out,0,sizeof(out));
   STARTTIMESHOT();
   fir_dsp(dataIn, dataOut, size, gauss, 23);
   STOPTIME("fir_dsp", 1);
   STARTTIMESHOT();
   for(int i=0; i<N; i++) fir_dsp(dataIn, dataOut, size, gauss, 23);
   STOPTIME("fir_dsp time",N);
 
   std::memset(out,0,sizeof(out));
   STARTTIMESHOT();
   fir_1_16(dataIn, (q15*)gauss, 24, dataOut, size);
   STOPTIME("fir_1_16", 1);
   // TRACEDATA("out",dataOut,50);
   // TRACEDATA("tail",&dataOut[size-10],34);
   STARTTIMESHOT();
   for(int i=0; i<N; i++) fir_1_16(dataIn, (q15*)gauss, 24, dataOut, size);
   STOPTIME("fir_1_16 time",N);

   std::memset(out,0,sizeof(out));
   STARTTIMESHOT();
   fir_16_16_q15(dataIn, (q15*)gauss, 24, dataOut, size);
   STOPTIME("fir_16_16_q15", 1);
   // TRACEDATA("out",dataOut,50);
   // TRACEDATA("tail",&dataOut[size-10],34);
   STARTTIMESHOT();
   for(int i=0; i<N; i++) fir_16_16_q15(dataIn, (q15*)gauss, 24, dataOut, size);
   STOPTIME("fir_16_16_q15 time",N);
 
   uint32_t mem2=esp_get_free_heap_size();
   if(mem1 != mem2)
   {
      TRACE("memory leak",mem1-mem2,false);
      TRACE("start",mem1,false);
      TRACE("stop",mem2,false);
      TEST_FAIL_MESSAGE("memory leak");
   }
}

#define toDegree(f) ((f*180)/M_PI)
TEST_CASE("complex", "[math][mathlib_s3]")
{
   uint32_t mem1=esp_get_free_heap_size();

   complex_q15 x = {toQ15(0.5),toQ15(0.5)};
   complex_q15 y = {toQ15(0.5),toQ15(0.5)};
   complex_q15 z = cmul_q15(x,y);
   std::printf("%d+i%d\n",z.re,z.im);
   y.im=-y.im;
   STARTTIMESHOT();
   z = cmul_q15(x,y);
   STOPTIME("cmul_q15", 1);
   //std::printf("%d+i%d\n",z.re,z.im);
   STARTTIMESHOT();
   for(int i=0; i<N; i++) z = cmul_q15(x,y);
   STOPTIME("cmul_q15", N);

   for(int i=0;i<10;i++)
   {
      sincos_q15(toQ15Angle((M_PI/4)*i), &in[i].im, &in[i].re);
   }
   std::memset(out,0,sizeof(out));
   complex_q15* cout = (complex_q15*)out;
   STARTTIMESHOT();
   cmul10_q15(in, &x, (complex_q15*)out);
   STOPTIME("cmul10_q15", 1);
   // TRACEDATA("in",(q15*)in,22);
   // TRACEDATA("out",out,22);
   // for(int i=0;i<10;i++)
   // {
   //    q15 angle = arg(in[i]);
   //    float a = toDegree(toFloatAngle(angle));
   //    std::printf(" %f",a);
   // }
   // std::printf("\n");
   // for(int i=0;i<10;i++)
   // {
   //    q15 angle = arg(cout[i]);
   //    float a = toDegree(toFloatAngle(angle));
   //    std::printf(" %f",a);
   // }
   // std::printf("\n");
   STARTTIMESHOT();
   for(int i=0; i<N; i++) cmul10_q15(in, &x, (complex_q15*)out);
   STOPTIME("cmul10_q15", N);
   std::printf("\n");

   uint32_t mem2=esp_get_free_heap_size();
   if(mem1 != mem2)
   {
      TRACE("memory leak",mem1-mem2,false);
      TRACE("start",mem1,false);
      TRACE("stop",mem2,false);
      TEST_FAIL_MESSAGE("memory leak");
   }
}
