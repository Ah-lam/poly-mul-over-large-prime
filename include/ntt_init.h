#pragma once

#include <stdlib.h>

#include "defs.h"
#include "mont_mul.h"

EXTERNC_BEGIN

typedef struct aligned128_ptr_s {
  void *    base;
  __uint128_t *ptr;
} aligned128_ptr_t;

typedef struct ntt_pre_table_s {
  uint64_t        m;        
  __uint128_t     w;
  __uint128_t w_inv; 
  mul_op256_t n_inv; 
  uint64_t        n;
  mul_op256_t    w1;

  aligned128_ptr_t w_powers;
  aligned128_ptr_t w_powers_con;
  aligned128_ptr_t w_inv_powers;
  aligned128_ptr_t w_inv_powers_con;

} ntt_pre_table_t;

int init_ntt_table();
void destroy_ntt_table();

const ntt_pre_table_t* get_ntt_table();

int table_len();

EXTERNC_END