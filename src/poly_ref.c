#include <string.h>
#include "poly_ref.h"


__uint128_t evaluate_poly(const __uint128_t co[], uint32_t len, __uint128_t x) {
  
  __uint128_t sum = co[0];
  __uint128_t mul = enter_mont(1);
  __uint128_t xm = x;

  for (uint32_t i = 1; i < len; i++) {
    mul = mul_mont(xm, mul);

    __uint128_t tmp = mul_mont(mul, co[i]);
    sum += tmp; 
  }
  return sum;
}


void polynomial_mul(const ntt_pre_table_t *t, const __uint128_t a[], const __uint128_t b[], __uint128_t c[]) {

  __uint128_t ar[t->n];
  __uint128_t br[t->n];

  size_t len = t->n >> 1;
  memcpy(ar, a, len * sizeof(__uint128_t));
  memcpy(br, b, len * sizeof(__uint128_t));

  memset(ar + len, 0, len * sizeof(__uint128_t));
  memset(br + len, 0, len * sizeof(__uint128_t));
  memset(c, 0, t->n * sizeof(__uint128_t));

  fwd_ntt(ar, t->n, t->w_powers.ptr, t->w_powers_con.ptr);

  fwd_ntt(br, t->n, t->w_powers.ptr, t->w_powers_con.ptr);

  for (uint64_t i = 0; i < t->n; i++) {
    c[i] = mul_mont(ar[i], br[i]);
  }

  inv_ntt(c, t->n, t->m, t->n_inv, t->w_inv_powers.ptr,
                     t->w_inv_powers_con.ptr, t->w1);

}

void polynomial_mul_normal(const ntt_pre_table_t *t, const __uint128_t a[], const __uint128_t b[], __uint128_t c[]) {
  // __uint128_t ar[t->n];
  // __uint128_t br[t->n];

  size_t len = t->n >> 1;
  // memcpy(ar, a, len * sizeof(__uint128_t));
  // memcpy(br, b, len * sizeof(__uint128_t));

  // memset(ar + len, 0, len * sizeof(__uint128_t));
  // memset(br + len, 0, len * sizeof(__uint128_t));
  memset(c, 0, t->n * sizeof(__uint128_t));
  
  for (size_t i = 0; i < len; i++) {
    for (size_t j = 0; j < len; j++) {
       c[i + j] += mul_mont(a[i], b[j]);
    }
  }
}


void poly_factory(poly_t p[], uint32_t poly_nums, uint32_t coff_len, __uint128_t buf[]) {

  uint32_t offset = 0;
  for (size_t i = 0; i < poly_nums; i++) {
    p[i].coef = &buf[offset];
    p[i].len = coff_len;
    offset += coff_len;
  }
}

void poly_clone(poly_t dst[], const poly_t src[], uint32_t poly_nums) {
  for (size_t i = 0; i < poly_nums; i++) {
    dst[i].coef = src[i].coef;
    dst[i].len = src[i].len;
  }
}


void poly_copy(poly_t dst[], const poly_t src[], uint32_t poly_nums, __uint128_t buf[]) {

  uint32_t coff_num = src[0].len;
  memset(buf, 0, sizeof(__uint128_t) * coff_num);
  uint32_t offset = 0;

  for (size_t i = 0; i < poly_nums; i++) {
    dst[i].coef = &buf[offset];
    memcpy(dst[i].coef, src[i].coef, sizeof(__uint128_t) * src[i].len);
    dst[i].len = src[i].len;
    offset += coff_num << 1;
  }
}

int poly_mul_twobytwo(const ntt_pre_table_t *t, const poly_t in[], uint32_t poly_nums, poly_t out[], polynomial_mul_function multiply) {

  size_t len = poly_nums >> 1;
  for (size_t i = 0; i < len; i ++) {
    multiply(t, in[i].coef, in[i + len].coef, out[i].coef);
  }
  return SUCCESS;
}

uint32_t bit_num(uint32_t n) {
    uint32_t num_bits = 0;
    while (n > 0) {
        n >>= 1;
        num_bits++;
    }
    return num_bits;
}

int poly_mul_continuous(const ntt_pre_table_t *t, const poly_t ins[], uint32_t poly_nums, poly_t *os) {

  __uint128_t ping[PINGPONG_BUFFER_SIZE];
  __uint128_t pong[PINGPONG_BUFFER_SIZE];
  poly_t in[poly_nums];
  poly_t out[poly_nums];
  memset(ping, 0, sizeof(__uint128_t) * PINGPONG_BUFFER_SIZE);
  memset(pong, 0, sizeof(__uint128_t) * PINGPONG_BUFFER_SIZE);

  poly_copy(in, ins, poly_nums, ping);

  uint32_t coff_nums = in[0].len;
  uint32_t next_coffs = 0;
  for (size_t num = poly_nums; num > 2; num = num >> 1) {
    next_coffs = coff_nums << 1;

    uint32_t next_poly_nums = num >> 1;
    poly_factory(out, next_poly_nums, next_coffs, pong);

    if (coff_nums <= 16) {
      ntt_pre_table_t st;
      st.n = next_coffs;

      poly_mul_twobytwo(&st, in, num, out, polynomial_mul_normal); 
    } else {
      uint32_t num_bits = bit_num(coff_nums);
      poly_mul_twobytwo(&t[num_bits - 4], in, num, out, polynomial_mul); 
    }

    poly_copy(in, out, next_poly_nums, ping); 
    coff_nums = next_coffs;
  }

  // last time:
  poly_factory(out, 1, coff_nums << 1, pong);
  if (coff_nums <= 16) {
    ntt_pre_table_t st;
    st.n = coff_nums << 1;
    poly_mul_twobytwo(&st, in, 2, out, polynomial_mul_normal); 
  } else {
    uint32_t num_bits = bit_num(coff_nums);
    poly_mul_twobytwo(&t[num_bits - 4], in, 2, out, polynomial_mul); 
  }

  memcpy(os->coef, out[0].coef, sizeof(__uint128_t) * out[0].len);
  os->len = out[0].len;

  return SUCCESS;
}
