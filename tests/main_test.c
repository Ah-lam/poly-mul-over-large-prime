#include "test.h"

#define NUM_TEST 500L

__uint128_t eval_poly(U128 co[], uint32_t len, U128 x) {
  
  __uint128_t sum = enter_mont(convert64to128(co[0]));
  __uint128_t mul = enter_mont(1);
  __uint128_t xm = enter_mont(convert64to128(x));

  for (uint32_t i = 1; i < len; i++) {
    __uint128_t com = enter_mont(convert64to128(co[i]));
    mul = mul_mont(xm, mul);
    __uint128_t tmp = mul_mont(mul, com);
    sum += tmp; 
  }
  return back_from_mont(sum);
}

void polynomial_mul_mont(const ntt_pre_table_t *t, const __uint128_t a[], const __uint128_t b[], __uint128_t c[]) {

  __uint128_t ar[t->n];
  __uint128_t br[t->n];

  memset(ar, 0, t->n * sizeof(__uint128_t));
  memset(br, 0, t->n * sizeof(__uint128_t));
  memset(c, 0, t->n * sizeof(__uint128_t));

  for (uint64_t i = 0; i < t->n / 2; i++) { 
    ar[i] = enter_mont(a[i]); 
    br[i] = enter_mont(b[i]);
  }

  fwd_ntt(ar, t->n, t->w_powers.ptr, t->w_powers_con.ptr);

  fwd_ntt(br, t->n, t->w_powers.ptr, t->w_powers_con.ptr);

  for (uint64_t i = 0; i < t->n; i++) {
    c[i] = mul_mont(ar[i], br[i]);
  }

  inv_ntt(c, t->n, t->m, t->n_inv, t->w_inv_powers.ptr,
                     t->w_inv_powers_con.ptr, t->w1);

  for (uint64_t i = 0; i < t->n; i++) {
    c[i] = back_from_mont(c[i]);
  }  

}

void polynomial_mul_standard(const ntt_pre_table_t *t, const __uint128_t a[], const __uint128_t b[], __uint128_t c[], __uint128_t q) {

  size_t len = t->n / 2;
  for (size_t i = 0; i < len; i++) {
    for (size_t j = 0; j < len; j++) {
       c[i + j] += mul_reduce(a[i], b[j]);
    }
  }

  for (uint64_t i = 0; i < t->n; i++) {
    c[i] = c[i] % q;
  }
}

int test_ntt(const ntt_pre_table_t *t) {

  __uint128_t q = convert64to128(Q);

  __uint128_t a_orig[t->n];
  random_buf(a_orig, t->n, q);
  __uint128_t a[t->n];

  memcpy(a, a_orig, sizeof(a));

  for (uint64_t i = 0; i < t->n; i++) {
    a[i] = enter_mont(a[i]);
  }

  fwd_ntt(a, t->n, t->w_powers.ptr, t->w_powers_con.ptr);
  inv_ntt(a, t->n, t->m, t->n_inv, t->w_inv_powers.ptr,
                     t->w_inv_powers_con.ptr, t->w1);

  for (uint64_t i = 0; i < t->n; i++) {
    a[i] = back_from_mont(a[i]);
  }  

  GUARD_MSG(memcmp(a_orig, a, sizeof(a)), "Bad results after ntt inv\n");
  return SUCCESS;
}

int test_polynomial_mul(const ntt_pre_table_t *t) {

  __uint128_t a[t->n];
  __uint128_t b[t->n];
  __uint128_t c1[t->n];
  __uint128_t c2[t->n];

  const __uint128_t q128 = convert64to128(Q);
  uint64_t len = t->n / 2;

  // Method 1: Using NTT calculation
  memset(a, 0, sizeof(a));
  memset(b, 0, sizeof(b));
  memset(c1, 0, sizeof(c1));

  random_buf(a, len, q128);       
  random_buf(b, len, q128);           

  polynomial_mul_mont(t, a, b, c1);

  // Method 2: Using standard calculation
  memset(c2, 0, sizeof(c2));

  polynomial_mul_standard(t, a, b, c2, q128);
  for (uint64_t i = 0; i < t->n; i++) {
    if (c1[i] != c2[i]) {
       printf("Bad results with standard calculating verification, c1[%lu] = 0x%016lx%016lx, c2[%lu] = 0x%016lx%016lx, \n", i, (uint64_t)(c1[i] >> 64), (uint64_t)c1[i], i, (uint64_t)(c2[i] >> 64), (uint64_t)c2[i]);
       return ERROR; 
    }
  }

  return SUCCESS;
}

int test_poly_mul_eval() {

  // 1. Generate simulated input polynomials
  const __uint128_t q128 = convert64to128(Q);

  U128 ping[MAX_NTT_NUM], pong[MAX_NTT_NUM];

  // The form of polynomial is (a + bx)
  const uint32_t coff_nums = 2; 

  const uint32_t poly_nums[3] = {8, 512, MAX_POLY_NUMS}; 
  for (uint32_t i = 0; i < sizeof(poly_nums)/ sizeof(uint32_t); i++) {
    memset(ping, 0, sizeof(U128) * MAX_NTT_NUM);
    memset(pong, 0, sizeof(U128) * MAX_NTT_NUM);

    random_buf_u128(ping, poly_nums[i] * coff_nums, q128);  

    r_poly_t ins[poly_nums[i]], outp;  
    poly_factory_u128(ins, poly_nums[i], coff_nums, ping); 

    uint32_t out_coff_num = ins[0].len << (bit_num(poly_nums[i]) - 1);
    poly_factory_u128(&outp, 1, out_coff_num, pong);

    // 2. Test API
    U128 x = {0}, gamma;
    random_buf_u128(&x, 1, q128); 

    GUARD(poly_mul_eval(ins, poly_nums[i], &outp, x, gamma))

    // Calculate the evaluation of entire polynomial
    __uint128_t r1 = eval_poly(outp.coef, out_coff_num, x);

    // 3. Comparison verification
    __uint128_t r2 = 1;
    for (size_t k = 0; k < poly_nums[i]; k++) {
      __uint128_t tmp = eval_poly(ins[k].coef, ins[k].len, x);
      r2 = mul_reduce(tmp, r2);
    }

    if (r1 != r2) {
      printf("failure to test test_poly_mul_eval case!\n");
      return ERROR;
    }
  }

  return SUCCESS;
}

int main() {

  srand((unsigned)time(NULL));

  GUARD(init_ntt_table())

  const ntt_pre_table_t *t = get_ntt_table();

  for (int i = 0; i < table_len(); i++) {
    for (int k = 0; k < NUM_TEST; k++) {
      GUARD(test_ntt(&t[i]))
      GUARD(test_polynomial_mul(&t[i]))
    }
  }

  for (int k = 0; k < NUM_TEST; k++) {
      GUARD(test_poly_mul_eval())
  }

  destroy_ntt_table();

  printf("case pass!\n");
  return SUCCESS;
}
