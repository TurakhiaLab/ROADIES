#pragma once

#include <cuda_runtime.h>
#include <cmath>

#if defined(MLIPPER_USE_DOUBLE)
using fp_t = double;
using fp2_t = double2;
using fp4_t = double4;
static constexpr const char* FP_MODE_NAME = "double";
#else
using fp_t = float;
using fp2_t = float2;
using fp4_t = float4;
static constexpr const char* FP_MODE_NAME = "float";
#endif

#if defined(MLIPPER_USE_DOUBLE)
static constexpr fp_t FP_EPS = fp_t(1e-300);
#else
static constexpr fp_t FP_EPS = fp_t(1e-30f);
#endif

__host__ __device__ __forceinline__ fp4_t make_fp4(fp_t x, fp_t y, fp_t z, fp_t w) {
#if defined(MLIPPER_USE_DOUBLE)
    return make_double4(x, y, z, w);
#else
    return make_float4(x, y, z, w);
#endif
}

__host__ __device__ __forceinline__ fp_t fp_exp(fp_t x) {
#if defined(MLIPPER_USE_DOUBLE)
    return ::exp(x);
#else
    return ::expf(x);
#endif
}

__host__ __device__ __forceinline__ fp_t fp_expm1(fp_t x) {
#if defined(MLIPPER_USE_DOUBLE)
    return ::expm1(x);
#else
    return ::expm1f(x);
#endif
}

__host__ __device__ __forceinline__ fp_t fp_log(fp_t x) {
#if defined(MLIPPER_USE_DOUBLE)
    return ::log(x);
#else
    return ::logf(x);
#endif
}

__host__ __device__ __forceinline__ fp_t fp_ldexp(fp_t x, int exp) {
#if defined(MLIPPER_USE_DOUBLE)
    return ::ldexp(x, exp);
#else
    return ::ldexpf(x, exp);
#endif
}

__host__ __device__ __forceinline__ fp_t fp_fabs(fp_t x) {
#if defined(MLIPPER_USE_DOUBLE)
    return ::fabs(x);
#else
    return ::fabsf(x);
#endif
}

__host__ __device__ __forceinline__ fp_t fp_fma(fp_t a, fp_t b, fp_t c) {
#if defined(MLIPPER_USE_DOUBLE)
    return ::fma(a, b, c);
#else
    return ::fmaf(a, b, c);
#endif
}

__host__ __device__ __forceinline__ fp_t fp_fmax(fp_t a, fp_t b) {
#if defined(MLIPPER_USE_DOUBLE)
    return ::fmax(a, b);
#else
    return ::fmaxf(a, b);
#endif
}

__host__ __device__ __forceinline__ fp_t fp_fmin(fp_t a, fp_t b) {
#if defined(MLIPPER_USE_DOUBLE)
    return ::fmin(a, b);
#else
    return ::fminf(a, b);
#endif
}

__host__ __device__ __forceinline__ fp_t fp_dot4(const fp4_t& a, const fp4_t& b) {
    return fp_fma(a.x, b.x,
           fp_fma(a.y, b.y,
           fp_fma(a.z, b.z, a.w * b.w)));
}

__host__ __device__ __forceinline__ fp_t fp_hmax4(fp_t a, fp_t b, fp_t c, fp_t d) {
    return fp_fmax(fp_fmax(a, b), fp_fmax(c, d));
}

__host__ __device__ __forceinline__ void fp_scale_pow2(fp_t& x, int shift) {
#if defined(__CUDA_ARCH__)
#if defined(MLIPPER_USE_DOUBLE)
    long long bits = __double_as_longlong(x);
    bits += ((long long)shift << 52);
    x = __longlong_as_double(bits);
#else
    int bits = __float_as_int(x);
    bits += (shift << 23);
    x = __int_as_float(bits);
#endif
#else
    x = fp_ldexp(x, shift);
#endif
}
