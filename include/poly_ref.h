#pragma once

#include "defs.h"
#include "mont_mul.h"
#include "ntt_init.h"
#include "ntt_ref.h"

EXTERNC_BEGIN

typedef struct poly_s {
  __uint128_t *coef;
  uint32_t len; 
} poly_t;

__uint128_t evaluate_poly(const __uint128_t co[], uint32_t len, __uint128_t x);
void polynomial_mul(const ntt_pre_table_t *t, const __uint128_t a[], const __uint128_t b[], __uint128_t c[]);
void polynomial_mul_normal(const ntt_pre_table_t *t, const __uint128_t a[], const __uint128_t b[], __uint128_t c[]);

typedef void (*polynomial_mul_function)(const ntt_pre_table_t *, const __uint128_t[], const __uint128_t[], __uint128_t[]);

void poly_factory(poly_t p[], uint32_t poly_nums, uint32_t coff_len, __uint128_t buf[]);
int poly_mul_twobytwo(const ntt_pre_table_t *t, const poly_t in[], uint32_t poly_nums, poly_t out[], polynomial_mul_function multiply);

uint32_t bit_num(uint32_t n);
int poly_mul_continuous(const ntt_pre_table_t *t, const poly_t ins[], uint32_t poly_nums, poly_t *os);

EXTERNC_END
