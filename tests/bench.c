#include "test.h"


static inline uint64_t cal_time(struct timespec time_start, struct timespec time_end) {
  return (time_end.tv_sec - time_start.tv_sec) * 1000000000 + (time_end.tv_nsec - time_start.tv_nsec);
}

#  define MEASURENTT(x, y) \
  const int total = 100000;                                                  \
  struct timespec time_start = { 0, 0 }, time_end = { 0, 0 };                     \
  clock_gettime(CLOCK_REALTIME, &time_start);                                     \
  for (uint64_t i = 0; i < total; i++) {                                          \
    x;                                                                            \
  }                                                                               \
  clock_gettime(CLOCK_REALTIME, &time_end);                                       \
  printf("%s: running time is %lu ns per opeartion.\n", y, cal_time(time_start, time_end) / total); 


int bench_poly_mul_eval() {

  const __uint128_t q128 = convert64to128(Q);

  U128 ping[MAX_NTT_NUM], pong[MAX_NTT_NUM];

  // The form of polynomial is (a + bx)
  const uint32_t coff_nums = 2; 
  const uint32_t poly_nums = MAX_POLY_NUMS; 

  memset(ping, 0, sizeof(U128) * MAX_NTT_NUM);
  memset(pong, 0, sizeof(U128) * MAX_NTT_NUM);

  random_buf_u128(ping, poly_nums * coff_nums, q128);  

  r_poly_t ins[poly_nums], outp;  
  poly_factory_u128(ins, poly_nums, coff_nums, ping); 

  uint32_t out_coff_num = ins[0].len << (bit_num(poly_nums) - 1);
  poly_factory_u128(&outp, 1, out_coff_num, pong);

  U128 x = {0}, gamma;
  random_buf_u128(&x, 1, q128); 

  /////////////////////////measurement/////////////////////
  MEASURENTT(poly_mul_eval(ins, poly_nums, &outp, x, gamma), "bench poly_mul_eval");

  return SUCCESS;          
}


int main() {

  srand((unsigned)time(NULL));

  GUARD(init_ntt_table())

  bench_poly_mul_eval();

  destroy_ntt_table();

  return SUCCESS;
}