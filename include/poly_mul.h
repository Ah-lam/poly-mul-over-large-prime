#pragma once

#include "defs.h"
#include "mont_mul.h"
#include "poly_ref.h"

EXTERNC_BEGIN

typedef struct r_poly_s {
  U128 *coef;
  uint32_t len; 
} r_poly_t;

int poly_mul_eval(const r_poly_t ins[], uint32_t poly_nums, r_poly_t *outp, const U128 x, U128 gamma);

EXTERNC_END