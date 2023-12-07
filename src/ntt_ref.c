#include "ntt_ref.h"


static inline void fwd_butterfly(__uint128_t *X, __uint128_t *Y, const mul_op256_t w)
{
  const __uint128_t q8 = convert64to128(Q) << 3;
  const __uint128_t X1 = *X;
  const __uint128_t T  = mul_mont(w.op, (*Y));
  *X = X1 + T;
  *Y = X1 - T + q8;
}

static inline void bkw_butterfly(__uint128_t *X, __uint128_t *Y, const mul_op256_t w, const __uint128_t q)
{

  const __uint128_t qn = q << 1;
  const __uint128_t X1 = *X + *Y;
  const __uint128_t T  = *X - *Y + qn;
  *X = X1;
  *Y = mul_mont(w.op, T);
}

static inline void bkw_butterfly_final(__uint128_t *     X,
                                       __uint128_t *     Y,
                                       const mul_op256_t w,
                                       const mul_op256_t n_inv, 
                                       const __uint128_t q)
{

  const __uint128_t qn = q << 2;
  const __uint128_t X1 = *X + *Y;
  const __uint128_t T  = *X - *Y + qn;
  *X = mul_mont(n_inv.op, X1);
  *Y = mul_mont(w.op, T);

}


void fwd_ntt(__uint128_t       a[],
             const uint64_t    N,
             const __uint128_t w[],
             const __uint128_t w_con[])
{
  size_t t = N >> 1;

  for(size_t m = 1; m < N; m <<= 1, t >>= 1) {
    size_t k = 0;
    for(size_t i = 0; i < m; i++) {
      const mul_op256_t w1 = {w[m + i], w_con[m + i]};

      for(size_t j = k; j < k + t; j++) {
        fwd_butterfly(&a[j], &a[j + t], w1);
      }
      k = k + (2 * t);
    }
  }
}

void inv_ntt(__uint128_t       a[],
             const uint64_t    N,
             const uint64_t    logN,
             const mul_op256_t n_inv,
             const __uint128_t w[],
             const __uint128_t w_con[],
             const mul_op256_t w1)
{
  uint64_t t = 1;
  const __uint128_t qn = convert64to128(Q) << logN;

  for(size_t m = N >> 1; m > 1; m >>= 1, t <<= 1) {
    size_t k = 0;
    for(size_t i = 0; i < m; i++) {
      const mul_op256_t w11 = {w[m + i], w_con[m + i]};

      for(size_t j = k; j < k + t; j++) {
        bkw_butterfly(&a[j], &a[j + t], w11, qn);
      }
      k = k + (2 * t);
    }
  }
 
  for(size_t j = 0; j < t; j++) {
    bkw_butterfly_final(&a[j], &a[j + t], w1, n_inv, qn);
  }
}

