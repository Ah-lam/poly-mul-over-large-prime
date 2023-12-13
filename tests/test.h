#pragma once

#include <stdint.h>
#include <string.h>
#include <time.h>
#include "../include/mont_mul.h"
#include "../include/poly_ref.h"
#include "../include/poly_mul.h"

EXTERNC_BEGIN

static inline uint64_t random64() {
    uint32_t h = random();
    uint32_t l = random();
    uint64_t x =  ((uint64_t)h << 32) + l ;
    return x;
}

static inline void random_buf(__uint128_t *values, const uint64_t n, const __uint128_t q) {
  for(uint64_t i = 0; i < n; i++) {
    uint64_t l = random64(); 
    uint64_t h = random64(); 

    values[i] = (__uint128_t)l | ((__uint128_t)h << 64);
  }
  for(uint64_t i = 0; i < n; i++) {
    values[i] = values[i] % q;
  }
}

static inline void random_buf_u128(U128 values[], const uint64_t n, const __uint128_t q) {
  for(uint64_t i = 0; i < n; i++) {
    uint64_t l = random64(); 
    uint64_t h = random64(); 

    __uint128_t tmp = (__uint128_t)l | ((__uint128_t)h << 64);
    tmp = tmp % q;
    convert128to64(tmp, values[i]);
  }
}

static inline void poly_factory_u128(r_poly_t p[], uint32_t poly_nums, uint32_t coff_len, U128 buf[]) {

  uint32_t offset = 0;
  for (size_t i = 0; i < poly_nums; i++) {
    p[i].coef = &buf[offset];
    p[i].len = coff_len;
    offset += coff_len;
  }
}

EXTERNC_END
