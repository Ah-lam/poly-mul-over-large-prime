#pragma once

#include "defs.h"
#include "mont_mul.h"

EXTERNC_BEGIN


void fwd_ntt(__uint128_t   a[],
                             uint64_t       N,
                             const __uint128_t w[],
                             const __uint128_t w_con[]);


void inv_ntt(__uint128_t       a[],
                        const uint64_t N,
                        const uint64_t    logN,
                        const mul_op256_t n_inv,
                        const __uint128_t w[],
                        const __uint128_t w_con[],
                        const mul_op256_t w1);

EXTERNC_END