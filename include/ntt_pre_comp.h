#pragma once

#include <string.h>

#include "defs.h"
#include "mont_mul.h"

EXTERNC_BEGIN

static inline __uint128_t bit_rev_idx(__uint128_t idx, uint64_t width)
{
  __uint128_t ret = 0;
  while(width > 0) {
    width--;
    ret += ((idx & 1) << width);
    idx >>= 1;
  }

  return ret;
}

static inline void bit_rev(__uint128_t       w_powers[],
                           const __uint128_t w[],
                           const __uint128_t N,
                           const uint64_t width)
{
  for(size_t i = 0; i < N; i++) {
    w_powers[bit_rev_idx(i, width)] = w[i];
  }
}

static inline void calc_w(__uint128_t       w_powers_rev[],
                          const __uint128_t w,
                          const __uint128_t N,
                          const uint64_t width)
{
  __uint128_t w_powers[N];
  w_powers[0] = 1;
  for(size_t i = 1; i < N; i++) {
    w_powers[i] = mul_reduce(w_powers[i - 1], w);
  }
  
  bit_rev(w_powers_rev, w_powers, N, width);
}

static inline void calc_w_inv(__uint128_t       w_inv_rev[],
                              const __uint128_t w_inv,
                              const __uint128_t N,
                              const uint64_t width)
{
  __uint128_t w_inv_powers[N];
  w_inv_powers[0] = 1;
  for(size_t i = 1; i < N; i++) {
    w_inv_powers[i] = mul_reduce(w_inv_powers[i - 1], w_inv);
  }

  bit_rev(w_inv_rev, w_inv_powers, N, width);
}

static inline void calc_w_con(__uint128_t       w_con[],
                              const __uint128_t w[],
                              const __uint128_t N)
{
  for(size_t i = 0; i < N; i++) {
    U256 wtmp = {0, 0, w[i] & MASK, (w[i] >> 64) & MASK};
    w_con[i] = divide_256(wtmp);
  }
}

static inline void expand_w(__uint128_t       w_expanded[],
                            const __uint128_t w[],
                            const __uint128_t N,
                            const __uint128_t q)
{
  w_expanded[0] = w[0];
  w_expanded[1] = 0;
  w_expanded[2] = w[1];
  w_expanded[3] = 0;
  for(size_t i = 4; i < 2 * N; i += 2) {
    w_expanded[i] = w[i / 2];

    if(i % 4 == 0) {
      const __uint128_t t = w_expanded[i / 2];
      w_expanded[i + 1] = mul_reduce(t, w[i / 2]); 
    } else {
      const __uint128_t t = w_expanded[(i - 2) / 2];
      w_expanded[i + 1] = q - mul_reduce(t, w[i / 2]); 
    }
  }
}


EXTERNC_END