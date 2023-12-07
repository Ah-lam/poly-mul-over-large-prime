#pragma once

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "domain.h"

EXTERNC_BEGIN

static const U128 Q = {1, (271 * (1 << 20))};

__uint128_t mul_reduce(__uint128_t a, __uint128_t b);

__uint128_t enter_mont(__uint128_t x);
__uint128_t mul_mont(__uint128_t a, __uint128_t b);
__uint128_t back_from_mont(__uint128_t x);

__uint128_t divide_256(U256 dividend);

void enter_mont_u128(U128 x, U128 r);
void back_from_mont_u128(U128 x, U128 r);

void accumu_mul_mont(U128 x[], uint32_t len, U128 r);

EXTERNC_END