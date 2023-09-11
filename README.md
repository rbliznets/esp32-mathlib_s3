# Библитотека алгориnмов (Q15) оптимизированных под esp32-s3
Для добавления в проект в папке компонентов из командной строки запустить:    

    git submodule add https://github.com/rbliznets/esp32-mathlib_s3 mathlib_s3
## БПФ
Реализация Radix2 с автонормированием и без. Время выполнения:
```
This is esp32s3 (240MHz) chip with 2 CPU core(s), WiFi/BLE, silicon revision 1, Minimum free heap size: 378624 bytes

**** FFT size 8192 *****
fft_radix2(1): 484usec
fft_radix2_scale(1): 904usec
revbin_permute(1): 260usec

**** FFT size 4096 *****
fft_radix2(1): 222usec
fft_radix2_scale(1): 376usec
revbin_permute(1): 130usec

**** FFT size 2048 *****
fft_radix2(1): 107usec
fft_radix2_scale(1): 183usec
revbin_permute(1): 66usec
```
## CORDIC
- atan2
- sinncos
```
sincos_q15(1000): 915nsec
```
## Complex
- arg
- magnitude
- cmul
```
**** Data size 1016 *****
arg_q15(1): 214usec
arg_q15(1000): 198usec
**** Data size 1 *****
float std::arg(1000): 1096nsec
arg_fr16(1000): 13usec
arg(1000): 709nsec
```
## Фильтр
```
**** Data size 1016 *****
dot_product_16_16(1000): 1765nsec
dot_product_1_16(1000): 2310nsec
dot_product_c(1000): 38usec

**** Data size 1000 *****
fir_1_16(1000): 196usec
fir_16_16_q15(1000): 36usec
```

