#include <string.h>
#include "poly_mul.h"

void r_poly_enter_mont(poly_t dst[], const r_poly_t src[], uint32_t poly_nums, __uint128_t buf[]) {

  uint32_t coff_num = src[0].len;
  memset(buf, 0, sizeof(__uint128_t) * coff_num);
  uint32_t offset = 0;

  for (size_t i = 0; i < poly_nums; i++) {
    dst[i].coef = &buf[offset];
    for (size_t k = 0; k < coff_num; k++) { 
      dst[i].coef[k] = enter_mont(convert64to128(src[i].coef[k])); 
    }
    dst[i].len = src[i].len;
    offset += coff_num;
  }
}

void poly_enter_mont(poly_t dst[], const poly_t src[], uint32_t poly_nums, __uint128_t buf[]) {

  uint32_t coff_num = src[0].len;
  memset(buf, 0, sizeof(__uint128_t) * coff_num);
  uint32_t offset = 0;

  for (size_t i = 0; i < poly_nums; i++) {
    dst[i].coef = &buf[offset];
    for (size_t k = 0; k < coff_num; k++) { 
      dst[i].coef[k] = enter_mont(src[i].coef[k]); 
    }
    dst[i].len = src[i].len;
    offset += coff_num;
  }
}

void poly_back_from_mont(const poly_t in[], r_poly_t out[], uint32_t poly_nums) {

  uint32_t coff_num = in[0].len;
  for (size_t i = 0; i < poly_nums; i++) {

    for (size_t k = 0; k < coff_num; k++) { 
      __uint128_t t = back_from_mont(in[i].coef[k]); 
      convert128to64(t, out[i].coef[k]);
    }
    out[i].len = in[i].len;
  }
}

static inline int is_power_of_two(uint32_t n) {
    if ((n & (n - 1)) == 0) {
        return 1;
    } else {
        return 0;
    }
}

int poly_mul_eval(const r_poly_t ins[], uint32_t poly_nums, r_poly_t *outp, const U128 x, U128 gamma) {

  if (poly_nums > MAX_POLY_NUMS || !is_power_of_two(poly_nums)) {
    return ERROR_ILLEGAL_PARAMETER;
  }

  const ntt_pre_table_t *t = get_ntt_table();
  if (t == NULL) {
    return ERROR_NTT_TABLE_EMPTY;
  }
  uint32_t out_coff_num = ins[0].len << (bit_num(poly_nums) - 1);
  if (out_coff_num > MAX_NTT_NUM) {
    return ERROR_OUT_OF_NTT_RANGE;
  }

  __uint128_t ping[PINGPONG_BUFFER_SIZE];
  poly_t in[poly_nums];
  memset(ping, 0, sizeof(__uint128_t) * PINGPONG_BUFFER_SIZE);

  poly_t outs;
  __uint128_t pong[PINGPONG_BUFFER_SIZE];
  memset(pong, 0, sizeof(__uint128_t) * PINGPONG_BUFFER_SIZE);
  
  poly_factory(&outs, 1, out_coff_num, pong);  

  r_poly_enter_mont(in, ins, poly_nums, ping);
  GUARD(poly_mul_continuous(t, in, poly_nums, &outs))

  __uint128_t xm = enter_mont(convert64to128(x));

  size_t half_len = outs.len >> 1;
  // only estimate polynomials of the third highest iterm
  __uint128_t r = evaluate_poly(outs.coef, half_len - 1, xm);

  __uint128_t gamma128 = back_from_mont(r);

  outs.len = half_len + 1;
  poly_back_from_mont(&outs, outp, 1);
  convert128to64(gamma128, gamma);

  return SUCCESS;
}