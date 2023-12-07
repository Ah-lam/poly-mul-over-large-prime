#include "mont_mul.h"
#include "divide_256.h"

// ====================================================
// constant values
//
// R: 0xDFFFFFFFFFFFF0E2B74312;
U128 R = {
    0xFFFFFFF0E2B74312,
    0xDFFFFF,
};

// R2 0x10646E78C139713E2EDB2345
U128 R2 = {
    0xC139713E2EDB2345,
    0x10646E78,
};

// ONE 0x01
U128 ONE = {
    0x1,
    0x0,
};


// ====================================================
// basic functions
//
// a*b
static inline void mul(uint64_t a, uint64_t b, U128 r) {
  __uint128_t a128 = (__uint128_t)(a);
  __uint128_t b128 = (__uint128_t)(b);
  __uint128_t r128 = a128 * b128;
  convert128to64(r128, r);
}

// a*b+c
static inline void mac(uint64_t a, uint64_t b, uint64_t c, U128 r) {
  __uint128_t a128 = (__uint128_t)(a);
  __uint128_t b128 = (__uint128_t)(b);
  __uint128_t c128 = (__uint128_t)(c);
  __uint128_t r128 = a128 * b128 + c128;
  convert128to64(r128, r);
}

// a*b+c+cin
static inline void macc(uint64_t a, uint64_t b, uint64_t c, uint64_t carry, U128 r) {
  __uint128_t a128 = (__uint128_t)(a);
  __uint128_t b128 = (__uint128_t)(b);
  __uint128_t c128 = (__uint128_t)(c);
  __uint128_t carry128 = (__uint128_t)(carry);
  __uint128_t r128 = a128 * b128 + c128 + carry128;
  convert128to64(r128, r);
}

// a+b
static inline void add(uint64_t a, uint64_t b, U128 r) {
  __uint128_t a128 = (__uint128_t)(a);
  __uint128_t b128 = (__uint128_t)(b);
  __uint128_t r128 = a128 + b128;
  convert128to64(r128, r);
}

// a+b+carry
static inline void adc(uint64_t a, uint64_t b, uint64_t carry, U128 r) {
  __uint128_t a128 = (__uint128_t)(a);
  __uint128_t b128 = (__uint128_t)(b);
  __uint128_t carry128 = (__uint128_t)(carry);
  __uint128_t r128 = a128 + b128 + carry128;
  convert128to64(r128, r);
}

// a-b
// borrow = -1 or 0
static inline void sub(uint64_t a, uint64_t b, U128 r) {
  __uint128_t a128 = (__uint128_t)(a);
  __uint128_t b128 = (__uint128_t)(b);
  __uint128_t r128 = a128 - b128;
  convert128to64(r128, r);
}

// a-b-borrow
// borrow = -1 or 0
static inline void sbb(uint64_t a, uint64_t b, uint64_t borrow, U128 r) {
  __uint128_t a128 = (__uint128_t)(a);
  __uint128_t b128 = (__uint128_t)(b);
  __uint128_t borrow128 = (__uint128_t)(borrow & 0x01);
  __uint128_t r128 = a128 - b128 - borrow128;
  convert128to64(r128, r);
}


// **********************************************************************
//
// u128 * u128
static inline void mul_128to256(U128 a, U128 b, U256 r) {
  U128 mul0, mul1, mul2, mul3;

  mul(a[0], b[0], mul0);                    // a0 * b0
  mac(a[0], b[1], mul0[1], mul1);           // a0 * b1
  mac(a[1], b[0], mul1[0], mul2);           // a1 * b0
  macc(a[1], b[1], mul1[1], mul2[1], mul3); // a1 * b1
  r[0] = mul0[0];
  r[1] = mul2[0];
  r[2] = mul3[0];
  r[3] = mul3[1];
}

// reduce
static inline void may_sub_q(uint64_t in0, uint64_t in1, uint64_t in2, U128 r) {
  U128 sub0, sub1, sub2;

  sub(in0, Q[0], sub0);
  sbb(in1, Q[1], sub0[1], sub1);
  sbb(in2, 0, sub1[1], sub2);

  if (sub2[1] == 0) {
    r[0] = sub0[0];
    r[1] = sub1[0]; 
  } else {
    r[0] = in0;
    r[1] = in1; 
  }
}

static inline void reduce(U256 in, U128 r) {

  U128 mul0;
  uint64_t k = 0 - in[0]; 
  mac(k, Q[0], in[0], mul0);
  macc(k, Q[1], in[1], mul0[1], mul0);

  U128 mul1;
  k = 0 - mul0[0];
  mac(k, Q[0], mul0[0], mul1);
  macc(k, Q[1], mul0[1], mul1[1], mul1);

  U128 add0, add1;
  add(in[2], mul1[0], add0);
  adc(in[3], mul1[1], add0[1], add1);
  may_sub_q(add0[0], add1[0], add1[1], r);
}

static inline void mul_mont_u128(U128 a, U128 b, U128 r) {
  U256 t;
  mul_128to256(a, b, t);
  reduce(t, r);
}

__uint128_t mul_mont(__uint128_t a, __uint128_t b) {

  U128 a128, b128, r128;
  convert128to64(a, a128);
  convert128to64(b, b128);

  U256 t;
  mul_128to256(a128, b128, t);
  reduce(t, r128);

  return convert64to128(r128);
}

__uint128_t enter_mont(__uint128_t x) {
  U128 x128 = {0}, r = {0};
  convert128to64(x, x128);
  mul_mont_u128(x128, R2, r);
  return convert64to128(r);
}

__uint128_t back_from_mont(__uint128_t x) {
  U128 x128 = {0}, r = {0};
  convert128to64(x, x128);
  mul_mont_u128(x128, ONE, r);
  return convert64to128(r);
}

__uint128_t mul_reduce(__uint128_t a, __uint128_t b) {
  U128 a128, b128;
  convert128to64(a, a128);
  convert128to64(b, b128);
  U128 r = {0};

  mul_mont_u128(a128, R2, a128);
  mul_mont_u128(b128, R2, b128);

  mul_mont_u128(a128, b128, r);

  mul_mont_u128(r, ONE, r);
  return convert64to128(r);
}


typedef struct pack_uint256_s {
  __uint128_t high;
  __uint128_t low;
} pack_uint256_t;

static inline pack_uint256_t sub_uint256(__uint128_t ah, __uint128_t al, __uint128_t bh, __uint128_t bl) {

  pack_uint256_t r;
  r.low = al - bl;
  r.high = ah - bh - (al < bl);
  return r;
}

static inline int cmp_uint256(__uint128_t ah, __uint128_t al, __uint128_t bh, __uint128_t bl) {

    if (ah != bh) {
        return (ah > bh) ? 1 : -1;
    } else if (al != bl) {
        return (al > bl) ? 1 : -1;
    } else {
        return 0;
    }
}
__uint128_t divide_256(U256 dividend) {

    __uint128_t quotient = 0;

    __uint128_t ldivisor = convert64to128(Q);
    __uint128_t endhigh = high128(dividend);
    __uint128_t endlow = low128(dividend);

    while (!(cmp_uint256(endhigh, endlow, 0, ldivisor) == -1)) {
      int shift = 0; 
      while (!(cmp_uint256(endhigh, endlow, high128(precomputed[shift]), low128(precomputed[shift])) == -1)) {
          shift++;
      }

      pack_uint256_t t = sub_uint256(endhigh, endlow, high128(precomputed[shift - 1]), low128(precomputed[shift - 1]));
      endhigh = t.high;
      endlow = t.low;
      quotient += low128(shiftone[shift - 1]);
    }
    return quotient;
}


void enter_mont_u128(U128 x, U128 r) { 
  mul_mont_u128(x, R2, r);
}

void back_from_mont_u128(U128 x, U128 r) {
  mul_mont_u128(x, ONE, r);
}

void accumu_mul_mont(U128 x[], uint32_t len, U128 r) { 
  
  r[0] = x[0][0];
  r[1] = x[0][1];
  for (uint32_t i = 1; i < len; i++) {
    mul_mont_u128(r, x[i], r);
  }
}