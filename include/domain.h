#pragma once

#include "defs.h"

EXTERNC_BEGIN

typedef struct mul_op256_s {
  __uint128_t op;
  __uint128_t con;   
} mul_op256_t;

typedef uint64_t U256[4] ALIGN(16); 
typedef uint64_t U128[2] ALIGN(16); 

// convert functions
//
// 128 -> (64, 64)
static inline void convert128to64(__uint128_t in, U128 r) {
  r[0] = in & MASK;
  r[1] = (in >> 64) & MASK;
  return;
}

// (64, 64)->128
static inline __uint128_t convert64to128(const U128 in) {
  return (__uint128_t)(in[0]) |
         ((__uint128_t)in[1] << 64);
}

static inline __uint128_t high128(U256 in) {
  return (__uint128_t)(in[2]) |
         ((__uint128_t)in[3] << 64);
}

static inline __uint128_t low128(U256 in) {
  return (__uint128_t)(in[0]) |
         ((__uint128_t)in[1] << 64);
}

EXTERNC_END