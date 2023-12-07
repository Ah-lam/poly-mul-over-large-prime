#include "ntt_init.h"
#include "ntt_pre_comp.h"

static ntt_pre_table_t *ntt_table = NULL;

typedef struct param_s {
  uint64_t m;    
  U128 w;
  U128 w_inv; // w^(-1) mod q
  U128 n_inv; // 2^(-m) mod q
  U128 n_invcon;
  U128 w1op;
  U128 w1con;
} param_t;

static const param_t params[] = {
  {.m = 4, .w = {0x0ef17d87cdae6d4a, 0x0000000000f046ff}, .w_inv = {0x1013ffb049a5b443, 0x0000000003c52f9a}, .n_inv = {0x0000000000000001, 0x000000000fe10000}, 
    .n_invcon = {0x00000000f1d48bce, 0xf000000000000000}, .w1op = {0x76558ea76c7ffdd1, 0x3079740}, .w1con = {0xe953f1e0f4fbe8d8, 0x2dca95482ad65b1f} },
  {.m = 5, .w = {0x0985e713df981e14, 0x00000000002d523a}, .w_inv = {0x9151196c54328d8e, 0x000000000d11b1e1}, .n_inv = {0x0000000000000001, 0x0000000010688000}, 
    .n_invcon = {0x0000000078ea45e7, 0xf800000000000000}, .w1op = {0xc4d538ac49c00118, 0x6f4345f}, .w1con = {0xb56070f85820b93, 0x691ab55bea94d270} },
  {.m = 6, .w = {0xbfe92ac68d636833, 0x0000000000714458}, .w_inv = {0x782895fc4dca0185, 0x000000000278a57d}, .n_inv = {0x0000000000000001, 0x0000000010ac4000}, 
    .n_invcon = {0x000000003c7522f3, 0xfc00000000000000}, .w1op = {0xe26a9c5624e0008c, 0x37a1a2f}, .w1con = {0x5ab0387c2c105c9, 0x348d5aadf54a6938} },
  {.m = 7, .w = {0x17d938e5696e1b6c, 0x00000000003b3fb1}, .w_inv = {0xbe0e9d812c98fcb1, 0x0000000007dea572}, .n_inv = {0x0000000000000001, 0x0000000010ce2000}, 
    .n_invcon = {0x000000001e3a9179, 0xfe00000000000000}, .w1op = {0xf1354e2b12700046, 0x1bd0d17}, .w1con = {0x2d581c3e16082e4, 0x1a46ad56faa5349c} },
  {.m = 8, .w = {0xd8635cfdc8b91602, 0x000000000029784a}, .w_inv = {0xa19d11437b25072f, 0x0000000001f7329e}, .n_inv = {0x0000000000000001, 0x0000000010df1000}, 
    .n_invcon = {0x000000000f1d48bc, 0xff00000000000000}, .w1op = {0xf89aa71589380023, 0xde868b}, .w1con = {0x16ac0e1f0b04172, 0xd2356ab7d529a4e}  },
  {.m = 9, .w = {0xbd9c3bdcb424a610, 0x000000000003b2e6}, .w_inv = {0xec04862b297fef73, 0x000000000e7fe0fd}, .n_inv = {0x0000000000000001, 0x0000000010e78800}, 
    .n_invcon = {0x00000000078ea45e, 0xff80000000000000}, .w1op = {0xfc4d538ac49c0012, 0x8e74345}, .w1con = {0xb56070f85820b9, 0x8691ab55bea94d27} },
  {.m = 10, .w = {0x9d66fc0839144f00, 0x0000000000012340}, .w_inv = {0xca1328d6acabc1a3, 0x0000000009831ca5}, .n_inv = {0x0000000000000001, 0x0000000010ebc400}, 
    .n_invcon = {0x0000000003c7522f, 0xffc0000000000000}, .w1op = {0x1d9563a9db1fff8, 0xc7c5e5d}, .w1con = {0x7fa54fc783d3efa3, 0xbcb72a5520ab596c} },
  {.m = 11, .w = {0x09331ad02a858ac3, 0x0000000000010659}, .w_inv = {0xa3d94e27557073bc, 0x000000000b450445}, .n_inv = {0x0000000000000001, 0x0000000010ede200}, 
    .n_invcon = {0x0000000001e3a917, 0xffe0000000000000}, .w1op = {0x7f1354e2b1270005, 0xab1d0d1}, .w1con = {0xc02d581c3e16082e, 0xa1a46ad56faa5349} },
};


static inline int allocate_aligned_array(aligned128_ptr_t *aptr, size_t qw_num) {
  size_t size_to_allocate = qw_num * sizeof(__uint128_t) + 64;
  if(NULL == ((aptr->base) = malloc(size_to_allocate))) {
    printf("Allocation error");
    return ERROR_MEMORY_ALLOCATION;
  }

  aptr->ptr = (__uint128_t *)(((uintptr_t)aptr->base & (~0x3fULL)) + 64);
  return SUCCESS;
}

static inline void free_aligned_array(aligned128_ptr_t *aptr)
{
  free(aptr->base);
  aptr->base = NULL;
  aptr->ptr  = NULL;
}

int init_ntt_table() {

  const size_t table_len = sizeof(params) / sizeof(param_t);

  // ntt_table
  if(NULL == (ntt_table = malloc(sizeof(ntt_pre_table_t) * table_len))) {
    printf("Allocation error");
    return ERROR_MEMORY_ALLOCATION;
  }

  for (size_t i = 0; i< table_len; i++) {
    ntt_table[i].m = params[i].m;
    ntt_table[i].n = 1UL << ntt_table[i].m;

    mul_op256_t n_inv;
    __uint128_t w = convert64to128(params[i].w);
    __uint128_t w_inv = convert64to128(params[i].w_inv);
    n_inv.op = convert64to128(params[i].n_inv);
    n_inv.con = convert64to128(params[i].n_invcon);

    // pre-compute tables 
    GUARD(allocate_aligned_array(&ntt_table[i].w_powers, ntt_table[i].n))
    calc_w(ntt_table[i].w_powers.ptr, w, ntt_table[i].n, ntt_table[i].m);

    GUARD(allocate_aligned_array(&ntt_table[i].w_powers_con, ntt_table[i].n))
    calc_w_con(ntt_table[i].w_powers_con.ptr, ntt_table[i].w_powers.ptr, ntt_table[i].n);

    GUARD(allocate_aligned_array(&ntt_table[i].w_inv_powers, ntt_table[i].n))
    calc_w_inv(ntt_table[i].w_inv_powers.ptr, w_inv, ntt_table[i].n, ntt_table[i].m);

    GUARD(allocate_aligned_array(&ntt_table[i].w_inv_powers_con, ntt_table[i].n))
    calc_w_con(ntt_table[i].w_inv_powers_con.ptr, ntt_table[i].w_inv_powers.ptr, ntt_table[i].n);

    // transfer all parameters to the montgomery domain
    ntt_table[i].w = enter_mont(w);
    ntt_table[i].w_inv = enter_mont(w_inv);
    ntt_table[i].n_inv.op = enter_mont(n_inv.op);
    ntt_table[i].n_inv.con = enter_mont(n_inv.con); 
    ntt_table[i].w1.op = enter_mont(convert64to128(params[i].w1op)); 
    ntt_table[i].w1.con = enter_mont(convert64to128(params[i].w1con)); 

    __uint128_t *w_powers_rev = ntt_table[i].w_powers.ptr;
    __uint128_t *w_con = ntt_table[i].w_powers_con.ptr;
    __uint128_t *w_inv_rev = ntt_table[i].w_inv_powers.ptr;
    __uint128_t *w_inv_con = ntt_table[i].w_inv_powers_con.ptr;
    for(size_t j = 0; j < ntt_table[i].n; j++) {
      w_powers_rev[j] = enter_mont(w_powers_rev[j]);
      w_con[j] = enter_mont(w_con[j]);
      w_inv_rev[j] = enter_mont(w_inv_rev[j]);
      w_inv_con[j] = enter_mont(w_inv_con[j]);
    }
  }

  return SUCCESS;
}

void destroy_ntt_table()
{

  for (size_t i = 0; i< sizeof(params) / sizeof(param_t); i++) {
    free_aligned_array(&ntt_table[i].w_powers);
    free_aligned_array(&ntt_table[i].w_powers_con);
    free_aligned_array(&ntt_table[i].w_inv_powers);
    free_aligned_array(&ntt_table[i].w_inv_powers_con);
  }
  free(ntt_table);
  ntt_table = NULL;
}


const ntt_pre_table_t* get_ntt_table() {
  return ntt_table;
}

int table_len() {
  return sizeof(params) / sizeof(param_t);
}
